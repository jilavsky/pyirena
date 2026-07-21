"""
    this contains needed hdf5 support for matilda
    used by saving blank BL_QRS data.
"""
import h5py
import os
import numpy as np
import six  #what is this for???
import logging


def readTextFile(path, filename, error_fraction=0.05):
    """
    Read text data files (.dat, .txt) with Q, Intensity, and optional Error columns.

    Files with only 2 columns (Q and I) are accepted. In that case the uncertainty
    is generated as  Error = Intensity * error_fraction  so that downstream tools
    always receive valid uncertainty data.

    Args:
        path: Directory path
        filename: File name
        error_fraction: Fraction of intensity used to generate uncertainty when no
                        error column is present in the file (default 0.05 = 5%).

    Returns:
        dict: Dictionary with 'Q', 'Intensity', 'Error', 'dQ' keys
              Returns None if file cannot be read
    """
    filepath = os.path.join(path, filename)

    try:
        # Try to read the file, detecting and skipping header rows
        # First, read all lines to find where numeric data starts
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Find the first line with numeric data
        skip_rows = 0
        for i, line in enumerate(lines):
            line = line.strip()
            if not line or line.startswith('#'):
                # Skip empty lines and comments
                skip_rows = i + 1
                continue
            # Try to parse first value as float
            try:
                parts = line.split()
                float(parts[0])
                # This is numeric data, start here
                skip_rows = i
                break
            except (ValueError, IndexError):
                # Not numeric, skip this line
                skip_rows = i + 1
                continue

        # Now read data with proper skip_rows
        data = np.loadtxt(filepath, comments='#', skiprows=skip_rows)

        if data.ndim == 1:
            # Single row, reshape
            data = data.reshape(1, -1)

        # Check number of columns
        if data.shape[1] < 2:
            logging.error(f"Text file {filename} must have at least 2 columns (Q, Intensity)")
            return None

        Q = data[:, 0]
        I = data[:, 1]

        # Check for error column; generate from intensity if absent
        if data.shape[1] >= 3:
            error = data[:, 2]
        else:
            error = I * error_fraction
            logging.info(
                f"No uncertainty column in {filename}; "
                f"generated as I × {error_fraction:.4f} (fractional uncertainty)"
            )

        # Check for dQ column
        dQ = None
        if data.shape[1] >= 4:
            dQ = data[:, 3]

        logging.info(f"Successfully read text file {filename} with {len(Q)} data points")
        return {
            'Q': Q,
            'Intensity': I,
            'Error': error,
            'dQ': dQ,
        }

    except Exception as e:
        logging.error(f"Error reading text file {filename}: {e}")
        return None


def readSimpleHDF5(path, filename):
    """
    Read simple HDF5 files with basic structure (e.g., /entry1/data1/Q and /entry1/data1/I).
    Fallback reader for files that don't follow full NXcanSAS standard.
    """
    Filepath = os.path.join(path, filename)
    try:
        with h5py.File(Filepath, 'r') as f:
            # Try common simple structures
            possible_paths = [
                ('entry1/data1', 'Q', 'I'),
                ('entry/data', 'Q', 'I'),
                ('data', 'Q', 'I'),
                ('', 'Q', 'I'),  # Root level
                ('', 'Q', 'Intensity'),
            ]

            for base, q_name, i_name in possible_paths:
                q_path = f"{base}/{q_name}".strip('/') if base else q_name
                i_path = f"{base}/{i_name}".strip('/') if base else i_name

                if q_path in f and i_path in f:
                    Q = f[q_path][()]
                    I = f[i_path][()]

                    # Try to find error data
                    error = None
                    dQ = None
                    for err_name in ['Idev', 'Error', 'error', 'I_error']:
                        err_path = f"{base}/{err_name}".strip('/') if base else err_name
                        if err_path in f:
                            error = f[err_path][()]
                            break

                    for dq_name in ['Qdev', 'dQ', 'Q_error']:
                        dq_path = f"{base}/{dq_name}".strip('/') if base else dq_name
                        if dq_path in f:
                            dQ = f[dq_path][()]
                            break

                    logging.info(f"Successfully read simple HDF5 file {filename} using path {base}")
                    return {
                        'Intensity': I,
                        'Q': Q,
                        'Error': error,
                        'dQ': dQ,
                    }

            logging.warning(f"No recognizable data structure found in {filename}")
            return None
    except Exception as e:
        logging.error(f"Error reading simple HDF5 file {filename}: {e}")
        return None


def _attr_first_str(value):
    """Normalise an HDF5 attribute value to a plain string.

    Handles scalars, bytes, and shape-(1,)/array attributes (NeXus often
    stores text attributes as 1-element arrays). For multi-element string
    attributes (e.g. ``I_axes = ['Q']`` or ``'Q,Q'``) the first token is
    returned, which is what 1D SAS data needs.
    """
    if value is None:
        return None
    # Unwrap numpy arrays / lists to their first element.
    if hasattr(value, 'ndim') and getattr(value, 'ndim', 0) >= 1:
        if value.size == 0:
            return None
        value = value.flat[0]
    elif isinstance(value, (list, tuple)):
        if not value:
            return None
        value = value[0]
    if isinstance(value, bytes):
        value = value.decode('utf-8', errors='replace')
    value = str(value)
    # Comma-separated axis lists -> first axis.
    if ',' in value:
        value = value.split(',')[0]
    return value.strip()


def _attr_all_str(value):
    """Normalise an HDF5 attribute to a list of tokens, splitting on commas.

    Unlike :func:`_attr_first_str` (which returns only the first token), this
    preserves every token — needed to detect the NXcanSAS slit-length
    declaration ``Q@resolutions = 'dQw,dQl'``, where ``dQl`` is the second
    token.
    """
    if value is None:
        return []
    if hasattr(value, 'ndim') and getattr(value, 'ndim', 0) >= 1:
        tokens = [t for t in np.asarray(value).ravel().tolist()]
    elif isinstance(value, (list, tuple)):
        tokens = list(value)
    else:
        tokens = [value]
    out = []
    for tok in tokens:
        if isinstance(tok, bytes):
            tok = tok.decode('utf-8', errors='replace')
        out.extend(part.strip() for part in str(tok).split(',') if part.strip())
    return out


def find_all_sasdata(hdf5_file):
    """Find every plottable SAS data group in an open HDF5 file.

    A SAS data group is any group with ``canSAS_class == 'SASdata'`` or any
    ``NXdata`` group carrying a ``signal`` attribute (the loose generic-NeXus
    convention used by some canSAS example files).

    Returns a list of dicts, one per dataset, in file order::

        {'path': '<hdf5 path to the NXdata group>',
         'entry': '<hdf5 path to the parent SASentry, or None>',
         'name': '<human-readable label>'}
    """
    found = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Group):
            return
        cls = _attr_first_str(obj.attrs.get('canSAS_class'))
        nx = _attr_first_str(obj.attrs.get('NX_class'))
        is_sasdata = (cls == 'SASdata') or (nx == 'NXdata' and 'signal' in obj.attrs)
        if not is_sasdata:
            return
        # Prefer an explicit canSAS_name, else the group's own basename.
        label = _attr_first_str(obj.attrs.get('canSAS_name')) or name.rsplit('/', 1)[-1]
        # Parent path is everything up to the last '/'.
        parent = name.rsplit('/', 1)[0] if '/' in name else None
        found.append({'path': name, 'entry': parent, 'name': label,
                      'canSAS_class': cls})

    hdf5_file.visititems(visit)
    return found


def _selectable_sasdata(f):
    """The list of SAS I(Q) curves a user would choose between.

    Restricts to groups explicitly marked ``canSAS_class == 'SASdata'`` so that
    auxiliary NXdata (raw fly scans, transmission spectra) are not offered as
    analysable curves. Falls back to every discovered data group for loose
    files that omit ``canSAS_class`` entirely.
    """
    all_data = find_all_sasdata(f)
    strict = [d for d in all_data if d.get('canSAS_class') == 'SASdata']
    return strict if strict else all_data


def _resolve_default_path(f):
    """Follow the NeXus ``@default`` attribute chain from the file root down to
    the canonical SAS data group.

    NXcanSAS uses ``@default`` at each level (root → NXentry → NXsubentry →
    NXdata) to mark the intended default plottable data. This walks that chain
    until it lands on a SAS data group, returning its path (or ``None`` if the
    chain does not resolve to one).
    """
    loc = ''  # '' == root
    for _ in range(10):  # guard against malformed/cyclic @default chains
        obj = f[loc] if loc else f
        if loc:
            cls = _attr_first_str(obj.attrs.get('canSAS_class'))
            nx = _attr_first_str(obj.attrs.get('NX_class'))
            if cls == 'SASdata' or (nx == 'NXdata' and 'signal' in obj.attrs):
                return loc
        nxt = _attr_first_str(obj.attrs.get('default'))
        if not nxt:
            break
        child = f"{loc}/{nxt}".strip('/')
        if child not in f:
            break
        loc = child
    return None


def _is_smr(d):
    """True when a dataset dict belongs to a Matilda slit-smeared (SMR) entry.

    Matilda always writes slit-smeared data under a group whose name ends with
    ``_SMR`` (e.g. ``entry/SampleName_SMR/sasdata``). Checking the HDF5 path
    is sufficient — the convention is stable because the USAXS pipeline is
    controlled by the same codebase.
    """
    path = d.get('path', '') or ''
    entry = d.get('entry', '') or ''
    return '_SMR/' in path or entry.endswith('_SMR')


def _filter_smr(datasets):
    """Return *datasets* with slit-smeared (SMR) entries removed.

    If all entries are SMR (shouldn't happen with well-formed files) the
    original list is returned unchanged so the caller always has something to
    load.
    """
    non_smr = [d for d in datasets if not _is_smr(d)]
    return non_smr if non_smr else datasets


def _ordered_sasdata(f, include_smr=False):
    """Selectable SAS curves with the ``@default`` dataset first.

    So ``[0]`` is always the file's canonical dataset (the one NXcanSAS marks
    ``@default`` — desmeared for Matilda USAXS files), while the rest preserve
    file order for the dataset picker.

    By default, slit-smeared copies (``_SMR`` entries written by Matilda) are
    stripped so that USAXS files presenting two copies of the same data never
    trigger the picker or the multi-dataset CLI warning.  Pass
    ``include_smr=True`` (used by the GUI "Slit smeared data" option) to keep
    them so the user can select the smeared curve.
    """
    selectable = _selectable_sasdata(f)
    datasets = selectable if include_smr else _filter_smr(selectable)
    default_path = _resolve_default_path(f)
    if default_path:
        for i, d in enumerate(datasets):
            if d['path'] == default_path and i != 0:
                datasets.insert(0, datasets.pop(i))
                break
    return datasets


def _smr_sibling_path(f, default_path):
    """Path to the ``_SMR`` slit-smeared sibling of a desmeared SASdata group.

    Matilda writes the desmeared curve under ``entry/<name>/sasdata`` and its
    slit-smeared twin under ``entry/<name>_SMR/sasdata``.  Given the desmeared
    (``@default``) path, return the smeared sibling's path, or ``None`` when the
    file has no slit-smeared entry.
    """
    if not default_path:
        # No resolved default: fall back to any SMR dataset present.
        smr = [d for d in _selectable_sasdata(f) if _is_smr(d)]
        return smr[0]['path'] if smr else None
    parts = default_path.rsplit('/', 1)
    if len(parts) == 2:
        entry, leaf = parts
        candidate = f"{entry}_SMR/{leaf}"
        if candidate in f:
            return candidate
    smr = [d for d in _selectable_sasdata(f) if _is_smr(d)]
    return smr[0]['path'] if smr else None


def _read_one_sasdata(f, data_path):
    """Read a single SASdata (NXdata) group into the standard ``Data`` dict.

    Robust to the two NXcanSAS axis conventions: the strict ``I_axes``
    attribute and the loose generic-NeXus ``axes`` attribute, falling back to
    a direct ``Q`` lookup. Likewise tolerates a missing ``signal`` attribute
    (falls back to ``I``) and missing ``uncertainties`` / ``resolutions``
    pointers (falls back to ``Idev`` / ``Qdev``).
    """
    grp = f[data_path]
    attributes = grp.attrs

    # Initialise all output variables so the return dict is always complete.
    intensity = None
    Int_attributes = {}
    units = None
    Kfactor = None
    OmegaFactor = None
    blankname = None
    thickness = None
    label = ""
    Q = None
    Q_attributes = {}
    Error = None
    Error_attributes = {}
    dQ = None
    dQ_attributes = {}

    # --- Intensity (signal) --------------------------------------------------
    signal_name = _attr_first_str(attributes.get('signal')) or 'I'
    i_path = f"{data_path}/{signal_name}"
    if i_path not in f and f"{data_path}/I" in f:
        i_path = f"{data_path}/I"          # last-resort fallback
    if i_path in f:
        dataset = f[i_path]
        intensity = dataset[()]
        Int_attributes = dataset.attrs
        units = Int_attributes.get('units')
        # Kfactor and other instrument-specific attributes are only present in
        # USAXS/calibrated data files; use .get() so plain NXcanSAS files read
        # without errors.
        Kfactor = Int_attributes.get("Kfactor")
        OmegaFactor = Int_attributes.get("OmegaFactor")
        blankname = Int_attributes.get("blankname")
        thickness = Int_attributes.get("thickness")
        label = Int_attributes.get("label", "")

    # --- Q axis: I_axes (strict NXcanSAS) | axes (generic NeXus) | 'Q' -------
    q_name = (_attr_first_str(attributes.get('I_axes'))
              or _attr_first_str(attributes.get('axes'))
              or 'Q')
    q_path = f"{data_path}/{q_name}"
    if q_path not in f and f"{data_path}/Q" in f:
        q_path = f"{data_path}/Q"          # last-resort fallback
    if q_path in f:
        dataset = f[q_path]
        Q = dataset[()]
        Q_attributes = dataset.attrs

    # --- Error (uncertainties pointer, else 'Idev') --------------------------
    uncertainties_key = _attr_first_str(Int_attributes.get('uncertainties')) or 'Idev'
    err_path = f"{data_path}/{uncertainties_key}"
    if err_path in f:
        dataset = f[err_path]
        Error = dataset[()]
        Error_attributes = dataset.attrs

    # --- Q resolution (resolutions pointer, else 'Qdev') ---------------------
    resolutions_key = _attr_first_str(Q_attributes.get('resolutions')) or 'Qdev'
    dq_path = f"{data_path}/{resolutions_key}"
    if dq_path in f:
        dataset = f[dq_path]
        dQ = dataset[()]
        dQ_attributes = dataset.attrs

    # --- Slit-smearing metadata (NXcanSAS slit-length resolution) ------------
    # Slit-smeared data declare Q@resolutions="dQw,dQl": the per-point width
    # dQw becomes dQ (above), and the scalar dQl is the slit (half-)length that
    # drives model smearing.  Presence of dQl -> the loaded curve IS slit
    # smeared; absence -> pinhole (slit_length 0.0).  This transparent,
    # attribute-driven detection matches how SasView decides to smear, so files
    # from any NXcanSAS pipeline (not just Matilda's ``_SMR`` naming) work.
    slit_length = 0.0
    resolution_tokens = _attr_all_str(Q_attributes.get('resolutions'))
    has_dql = 'dQl' in resolution_tokens or f"{data_path}/dQl" in f
    if has_dql and f"{data_path}/dQl" in f:
        try:
            slit_length = float(np.asarray(f[f"{data_path}/dQl"][()]).ravel()[0])
        except Exception:
            slit_length = 0.0
    is_slit_smeared = slit_length > 0

    return {
        'Intensity': intensity,
        'Q': Q,
        'dQ': dQ,
        'Error': Error,
        'units': units,
        'Int_attributes': Int_attributes,
        'Q_attributes': Q_attributes,
        'Error_Attributes': Error_attributes,
        'dQ_Attributes': dQ_attributes,
        "Kfactor": Kfactor,
        "OmegaFactor": OmegaFactor,
        "blankname": blankname,
        "thickness": thickness,
        'label': label,
        # Canonical slit-smearing fields carried end-to-end (see slit-smearing
        # plan §3.2): slit_length in 1/Å (0.0 => pinhole), is_slit_smeared bool.
        'slit_length': slit_length,
        'is_slit_smeared': is_slit_smeared,
    }


def list_nxcansas_datasets(path, filename, include_smr=False):
    """List every SAS data group in an NXcanSAS file (without loading arrays).

    Convenience wrapper around :func:`find_all_sasdata` that opens the file.
    Used by the GUI to offer a dataset picker and by the CLI to report how
    many datasets a multi-dataset file contains.

    Parameters
    ----------
    include_smr : bool
        When True, keep slit-smeared (``_SMR``) entries in the list so the GUI
        can offer them for selection.  Default False (desmeared only).

    Returns a list of ``{'path', 'entry', 'name'}`` dicts (see
    :func:`find_all_sasdata`).
    """
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r') as f:
        return _ordered_sasdata(f, include_smr=include_smr)


def file_has_smr_entry(path, filename):
    """True when the NXcanSAS file contains a slit-smeared (``_SMR``) dataset.

    Lets the GUI decide whether to show the "Slit smeared data" checkbox for a
    given file (only Matilda-style files carrying both a desmeared and a
    slit-smeared copy need it).
    """
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r') as f:
        return any(_is_smr(d) for d in _selectable_sasdata(f))


def readGenericNXcanSAS(path, filename, data_path=None, prefer_slit_smeared=False):
    """Read NXcanSAS data from a NeXus file.

    Only the minimum required to plot/analyse is needed: a Q/I(/Idev) triplet
    inside a SASdata group. The reader tolerates both the strict NXcanSAS
    ``I_axes`` attribute and the loose generic-NeXus ``axes`` attribute, and
    accepts SASentry groups written as either ``NXentry`` (Igor/Irena) or
    ``NXsubentry`` (pyirena).

    Args:
        path:      Directory containing the file.
        filename:  File name.
        data_path: Optional explicit HDF5 path to a specific SASdata group
                   (as returned by :func:`list_nxcansas_datasets`). When given,
                   exactly that dataset is read. When ``None`` and the file
                   holds more than one dataset, the FIRST one is read and a
                   warning is logged listing all available datasets — callers
                   that want user selection should call
                   :func:`list_nxcansas_datasets` first.

    Returns:
        dict with keys ``Intensity``, ``Q``, ``dQ``, ``Error`` (+ attribute
        sub-dicts), or ``None`` if no SAS data could be located. Falls back to
        :func:`readSimpleHDF5` for files with no recognisable NXcanSAS group.
    """
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r') as f:
        # If the caller already chose a dataset, just read it.
        if data_path is not None:
            if data_path in f:
                return _read_one_sasdata(f, data_path)
            logging.warning(f"Requested data path '{data_path}' not found in {filename}.")
            return None

        # Slit-smeared preference: load the _SMR sibling of the @default
        # dataset when one exists (used by the "Slit smeared data" option).
        if prefer_slit_smeared:
            default_path = _resolve_default_path(f)
            if default_path is None:
                ordered = _ordered_sasdata(f)
                default_path = ordered[0]['path'] if ordered else None
            smr_path = _smr_sibling_path(f, default_path)
            if smr_path is not None:
                return _read_one_sasdata(f, smr_path)
            logging.info(
                f"{filename}: prefer_slit_smeared=True but no slit-smeared "
                "(_SMR) entry found; loading the default (desmeared) dataset."
            )

        datasets = _ordered_sasdata(f)
        if not datasets:
            logging.warning(
                f"No NXcanSAS data groups found in the file {filename}. "
                "Trying simple HDF5 reader."
            )
            return readSimpleHDF5(path, filename)

        if len(datasets) > 1:
            names = ", ".join(d['name'] for d in datasets)
            logging.warning(
                f"{filename} contains {len(datasets)} data sets ({names}). "
                f"pyIrena loaded the first one ('{datasets[0]['name']}'). "
                "Use the GUI dataset picker or pass data_path= to choose another."
            )

        return _read_one_sasdata(f, datasets[0]['path'])

def saveNXcanSAS(Sample,path, filename):
    
    #read stuff from the data dictionary
    Intensity = Sample["CalibratedData"]["Intensity"]
    Q = Sample["CalibratedData"]["Q"]
    Error = Sample["CalibratedData"]["Error"]
    dQ = Sample["CalibratedData"]["dQ"]
    Kfactor = Sample["CalibratedData"]["Kfactor"] if "Kfactor" in Sample["CalibratedData"] else None
    OmegaFactor = Sample["CalibratedData"]["OmegaFactor"] if "OmegaFactor" in Sample["CalibratedData"] else None
    blankname = Sample["CalibratedData"]["blankname"] if "blankname" in Sample["CalibratedData"] else None
    thickness = Sample["CalibratedData"]["thickness"] if "thickness" in Sample["CalibratedData"] else None
    label = Sample["RawData"]["filename"]
    timeStamp = Sample["RawData"]["metadata"]["timeStamp"]
    samplename = Sample["RawData"]["sample"]["name"]
    if isinstance(samplename, bytes):
        samplename = samplename.decode('utf-8')


    if "SMR_Int" in Sample["CalibratedData"]:
        SMR_Int =Sample["CalibratedData"]["SMR_Int"]
        SMR_Error =Sample["CalibratedData"]["SMR_Error"]
        SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
        SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]
        slitLength=Sample["CalibratedData"]["slitLength"]
    else:
        SMR_Int = None
        SMR_Error = None
        SMR_Qvec = None
        SMR_dQ = None
        slitLength = None

    R_Int = Sample["reducedData"]["Intensity"]
    R_Qvec = Sample["reducedData"]["Q"]
    R_Error=Sample["reducedData"]["Error"]

    if "BlankData" in Sample:
        BL_R_Int = Sample["BlankData"]["Intensity"]
        BL_Q_vec = Sample["BlankData"]["Q"]
        BL_Error = Sample["BlankData"]["Error"]    
    else:
        BL_R_Int = None
        BL_Q_vec = None
        BL_Error = None
        
    #this is Desmeared USAXS data, SLitSmeared data and plot data, all at once.
    # create the HDF5 NeXus file with same structure as our raw data files have...
    Filepath = os.path.join(path, filename)
    logging.info(f"Saving NXcanSAS data to {Filepath}")
    with h5py.File(Filepath, "a") as f:
        # point to the default data to be plotted
        f.attrs['default']          = 'entry'   #our files have one entry input.
        # these are hopefully optional and useful. 
        f.attrs['file_name']        = filename
        f.attrs['file_time']        = timeStamp 
        f.attrs['instrument']       = '12IDE USAXS'
        f.attrs['creator']          = 'Matilda NeXus writer'
        f.attrs['Matilda_version']  = '1.0.0' # version 2025-07-06
        f.attrs['NeXus_version']    = '4.3.0' #2025-5-9 4.3.0 is rc, it is current. 
        f.attrs['HDF5_version']     = six.u(h5py.version.hdf5_version)
        f.attrs['h5py_version']     = six.u(h5py.version.version)

        # now create the NXentry group called entry if does not exist
        if 'entry' not in f:
            nxentry = f.create_group('entry')    
        
        nxentry = f['entry']
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['canSAS_class'] = 'SASentry'
        nxentry.attrs['default']  = samplename   #modify with the most reduced data.
        
        #add definition as NXsas - this is location of raw AND reduced data
        # Check if 'definition' dataset exists in the entry group and delete it if present
        if 'definition' in nxentry:
            del nxentry['definition']
        nxentry.create_dataset('definition', data='NXsas')
        # other groups should be here from RAW data, so ignore. 

        if Intensity is not None:
            logging.info(f"Wrote Desmeared NXcanSAS group for file {filename}. ")
            # create the NXsubentry group for Desmeared reduced data. 
            newDataPath = "entry/"+samplename
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]
            
            nxDataEntry = f.create_group(newDataPath)
            nxDataEntry.attrs['NX_class'] = 'NXsubentry'
            nxDataEntry.attrs['canSAS_class'] = 'SASentry'
            nxDataEntry.attrs['default'] = 'sasdata'
            nxDataEntry.attrs['title'] = samplename
            #add definition as NXcanSas
            nxDataEntry.create_dataset('definition', data='NXcanSAS')
            #add title as NXcanSas
            nxDataEntry.create_dataset('title', data=samplename)
            #add run (compulsory)
            nxDataEntry.create_dataset('run', data="run_identifier")

            # create the NXdata group for I(Q) for the avergaed data
            nxdata = nxDataEntry.create_group('sasdata')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'Q'      # X axis of default plot
            #nxdata.attrs['Q_indices'] = [1]    # TODO not sure what this means

            # Y axis data
            ds = nxdata.create_dataset('I', data=Intensity)
            ds.attrs['units'] = '1/cm'
            ds.attrs['uncertainties'] = 'Idev'
            ds.attrs['long_name'] = 'cm2/cm3'    # suggested X axis plot label
            ds.attrs['blankname'] = blankname
            ds.attrs['thickness'] = thickness
            ds.attrs['label'] = label
            ds.attrs['long_name'] = 'Intensity'    # suggested X axis plot label
            if Kfactor is not None:
                ds.attrs['Kfactor'] = Kfactor
            if OmegaFactor is not None:
                ds.attrs['OmegaFactor'] = OmegaFactor

            # X axis data
            ds = nxdata.create_dataset('Q', data=Q)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label
            ds.attrs['resolutions'] = 'Qdev'
        
            # d X axis data
            ds = nxdata.create_dataset('Qdev', data=dQ)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'   
            # dI axis data
            ds = nxdata.create_dataset('Idev', data=Error)
            ds.attrs['units'] = 'cm2/cm3'
            ds.attrs['long_name'] = 'Uncertainties'  

        if SMR_Int is not None:
            logging.info(f"Wrote SMR NXcanSAS group for file {filename}. ")
            # add the SMR data
            # create the NXsubentry group for Desmeared reduced data. 
            newDataPath = "entry/"+samplename+"_SMR"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            nxDataEntry.attrs['NX_class'] = 'NXsubentry'
            nxDataEntry.attrs['canSAS_class'] = 'SASentry'
            nxDataEntry.attrs['default'] = 'sasdata'
            nxDataEntry.attrs['title'] = samplename
            #add definition as NXcanSas
            nxDataEntry.create_dataset('definition', data='NXcanSAS')
            #add title as NXcanSas
            nxDataEntry.create_dataset('title', data=samplename)
            #add run (compulsory)
            nxDataEntry.create_dataset('run', data="run_identifier")

            # create the NXdata group for I(Q) for the avergaed data
            nxdata = nxDataEntry.create_group('sasdata')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'Q'      # X axis of default plot
            #nxdata.attrs['Q_indices'] = [1]    # TODO not sure what this means

            # Y axis data
            ds = nxdata.create_dataset('I', data=SMR_Int)
            ds.attrs['units'] = '1/cm'
            ds.attrs['uncertainties'] = 'Idev'
            ds.attrs['long_name'] = 'Intensity[cm2/cm3]'    # suggested X axis plot label
            ds.attrs['Kfactor'] = Kfactor
            ds.attrs['OmegaFactor'] = OmegaFactor
            ds.attrs['blankname'] = blankname
            ds.attrs['thickness'] = thickness
            ds.attrs['label'] = label

            # X axis data
            ds = nxdata.create_dataset('Q', data=SMR_Qvec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label
            ds.attrs['resolutions'] = 'dQw,dQl'
        
            # d X axis data
            ds = nxdata.create_dataset('dQw', data=SMR_dQ)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'dQw (A^-1)'           
            # slitlength
            ds = nxdata.create_dataset('dQl', data=slitLength)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'dQl (A^-1)'   
            # dI axis data
            ds = nxdata.create_dataset('Idev', data=SMR_Error)
            ds.attrs['units'] = 'cm2/cm3'
            ds.attrs['long_name'] = 'Uncertainties'  

        if R_Int is not None:
            logging.info(f"Wrote QRS group for file {filename}. ")
            newDataPath = "entry/"+"QRS_data"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            # R_Int axis data
            ds = nxDataEntry.create_dataset('Intensity', data=R_Int)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Intensity'    # suggested X axis plot label
            # R_Qvec axis data
            ds = nxDataEntry.create_dataset('Q', data=R_Qvec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q'    # suggested X axis plot label
            # R_Error axis data
            ds = nxDataEntry.create_dataset('Error', data=R_Error)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Error'    # suggested X axis plot label
            
        if BL_R_Int is not None:
            logging.info(f"Wrote Blank data group for file {filename}. ")
            newDataPath = "entry/"+"Blank_data"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            # R_Int axis data
            ds = nxDataEntry.create_dataset('Intensity', data=BL_R_Int)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Intensity'    # suggested X axis plot label
            ds.attrs['blankname']=blankname 
            # R_Qvec axis data
            ds = nxDataEntry.create_dataset('Q', data=BL_Q_vec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q'    # suggested X axis plot label
            # R_Error axis data
            ds = nxDataEntry.create_dataset('Error', data=BL_Error)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Error'    # suggested X axis plot label
 
 

    logging.info(f"Wrote NXcanSAS data to file: {filename}")

def readMyNXcanSAS(path, filename, isUSAXS = False):
    """
    Read My own data from NXcanSAS data in Nexus file.
    
    Parameters:
    path (str): The directory path where the file is located.
    filename (str): The name of the Nexus file to read.

    Returns:
    dict: A dictionary containing the read data.
    """    
    Sample = dict()
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r') as f:
        # Start at the root
        # Find the NXcanSAS entries 
        required_attributes = {'canSAS_class': 'SASentry', 'NX_class': 'NXsubentry'}
        required_items = {'definition': 'NXcanSAS'}
        SASentries =  find_matching_groups(f, required_attributes, required_items)
        logging.debug(f"Found {SASentries} entries in the file:{filename}")
              
        location = 'entry/QRS_data/'
        if location in f:
            Sample['reducedData'] = dict()
            dataset = _get_h5_value(f, location + "Intensity")
            if dataset is not None:
                Sample['reducedData']['Intensity'] = dataset
            dataset = _get_h5_value(f, location + "Q")
            if dataset is not None:
                Sample['reducedData']['Q'] = dataset
            dataset = _get_h5_value(f, location + "Error")
            if dataset is not None:
                Sample['reducedData']['Error'] = dataset

        location = 'entry/Blank_data/'
        if location in f:
            Sample['BlankData'] = dict()
            dataset = _get_h5_value(f, location + "Intensity")
            if dataset is not None:
                Sample['BlankData']['Intensity'] = dataset
            dataset = _get_h5_value(f, location + "Q")
            if dataset is not None:
                Sample['BlankData']['Q'] = dataset
            dataset = _get_h5_value(f, location + "Error")
            if dataset is not None:
                Sample['BlankData']['Error'] = dataset
            # BL_R_Int = Sample["BlankData"]["Intensity"]
            # BL_Q_vec = Sample["BlankData"]["Q"]
            # BL_Error = Sample["BlankData"]["Error"]    

        #location = 'entry/'+filename.split('.')[0]+'_SMR/'
        # location is the first of entries from SASentries which contains string _SMR
        location = next((entry + '/' for entry in SASentries if '_SMR' in entry), None)
        logging.debug(f"Found SMR entry at: {location}")
        if 'CalibratedData' not in Sample:
            Sample['CalibratedData'] = dict()
        
        if location is not None and location in f:
            isUSAXS = True      #have SMR data, assume USAXS setup

            dataset = _get_h5_value(f, location + "sasdata/I")
            if dataset is not None:
                Sample['CalibratedData']['SMR_Int'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Q")
            if dataset is not None:
                Sample['CalibratedData']['SMR_Qvec'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Idev")
            if dataset is not None:
                Sample['CalibratedData']['SMR_Error'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/dQw")
            if dataset is not None:
                Sample['CalibratedData']['SMR_dQ'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/dQl")
            if dataset is not None:
                Sample['CalibratedData']['slitLength'] = dataset
        else:
            Sample["CalibratedData"] ["SMR_Qvec"] = None,
            Sample["CalibratedData"] ["SMR_Int"] = None,
            Sample["CalibratedData"] ["SMR_Error"] = None,
            Sample["CalibratedData"] ["SMR_dQ"] = None,
            Sample["CalibratedData"] ["slitLength"] = None,

     
        location = next((entry + '/' for entry in SASentries if '_SMR' not in entry), None)
        logging.debug(f"Found NXcanSAS entry at: {location}")
        if 'RawData' not in Sample:          
            Sample['RawData'] = dict()
            Sample["RawData"]["sample"]=dict()

        if location is not None and location in f:
            dataset = _get_h5_value(f, location + "sasdata/I")
            if dataset is not None:
                Sample['CalibratedData']['Intensity'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Q")
            if dataset is not None:
                Sample['CalibratedData']['Q'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Idev")
            if dataset is not None:
                Sample['CalibratedData']['Error'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Qdev")
            if dataset is not None:
                Sample['CalibratedData']['dQ'] = dataset
            dataset = _get_h5_value(f, location + "title")
            if dataset is not None:
                Sample["RawData"]["sample"]["name"] = dataset

            location = location+'/sasdata/'
            if "I" in f[location]:
                attributes = f[location + "I"].attrs
                Sample['CalibratedData']['units'] = attributes['units']
                Sample['CalibratedData']['blankname'] = attributes["blankname"]
                Sample['CalibratedData']['thickness'] = attributes["thickness"]
                Sample["RawData"]["filename"] = attributes["label"]
                Sample['CalibratedData']['Kfactor'] = attributes["Kfactor"] if "Kfactor" in attributes else None
                Sample['CalibratedData']['OmegaFactor'] = attributes["OmegaFactor"] if "OmegaFactor" in attributes else None
        else:
            Sample["RawData"]["filename"] = filename
            Sample['CalibratedData']['Intensity'] = None
            Sample['CalibratedData']['Q'] = None
            Sample['CalibratedData']['Error'] = None
            Sample['CalibratedData']['Kfactor'] = None
            Sample['CalibratedData']['OmegaFactor'] = None
            Sample['CalibratedData']['blankname'] = None
            Sample['CalibratedData']['thickness'] = None
            Sample['CalibratedData']['units'] = None
            Sample['CalibratedData']['Error'] = None

        #and now we need to read the other groups, which are raw data... 
        if isUSAXS : 
            #metadata
            keys_to_keep = ['AR_center', 'ARenc_0', 'DCM_energy', 'DCM_theta', 'I0Gain','detector_distance',
                            'timeStamp','I0AmpGain',
                            'trans_pin_counts','trans_pin_gain','trans_pin_time','trans_I0_counts','trans_I0_gain',
                            'UPDsize', 'trans_I0_counts', 'trans_I0_gain', 'upd_bkg0', 'upd_bkg1','upd_bkg2','upd_bkg3',
                            'upd_bkgErr0','upd_bkgErr1','upd_bkgErr2','upd_bkgErr3','upd_bkgErr4','upd_bkg_err0',
                            'upd_bkg4','DDPCA300_gain0','DDPCA300_gain1','DDPCA300_gain2','DDPCA300_gain3','DDPCA300_gain4',
                            'SAD_mm', 'SDD_mm', 'thickness', 'title', 'useSBUSAXS',
                            'intervals', 'VToFFactor',
                            'upd_amp_change_mask_time0','upd_amp_change_mask_time1','upd_amp_change_mask_time2','upd_amp_change_mask_time3','upd_amp_change_mask_time4',
                        ]
            # Prefer the "classic" location, but fall back to the Bluesky one.
            metadata_group = None
            for md_path in (
                    "/entry/metadata",
                    "/entry/instrument/bluesky/metadata",
            ):
                if md_path in f:
                    metadata_group = f[md_path]
                    break
            if metadata_group is None:
                raise KeyError(
                    f"Could not find metadata group in {filename}. "
                    "Tried: /entry/metadata and /entry/instrument/bluesky/metadata"
                )

            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            # we need this key to be there also... Copy of the other one.
            # Ensure I0AmpGain exists; default to 1e6 if missing.
            metadata_dict["I0AmpGain"] = metadata_dict.get("I0AmpGain", 1e6)
            metadata_dict["I0Gain"] = metadata_dict["I0AmpGain"]
            #Instrument
            keys_to_keep = ['monochromator', 'energy', 'wavelength']
            instrument_group = f['/entry/instrument']
            instrument_dict = read_group_to_dict(instrument_group)
            instrument_dict = filter_nested_dict(instrument_dict, keys_to_keep)
            # sample
            sample_group = f['/entry/sample']
            sample_dict = read_group_to_dict(sample_group)

            Sample["RawData"]["metadata"] = metadata_dict
            Sample["RawData"]["instrument"] = instrument_dict
            Sample["RawData"]["sample"].update(sample_dict)
        else:       #this is SWAXS
            #metadata
            instrument_group = f['/entry/instrument']
            instrument_dict = read_group_to_dict(instrument_group)
            #occasionally this fails since 'data' does not exist. 
            # now, why this shoudl tno exists is mystery for me... 
            try:
                del instrument_dict['detector']['data']
            except KeyError:
                pass
            #metadata
            keys_to_keep = ['I000_cts', 'I00_cts', 'I00_gain', 'I0_cts', 'I0_cts_gated',
                            'TR_cts_gated','TR_cts','TR_gain','I0_Sample',
                            'I0_gain', 'I_scaling', 'Pin_TrI0', 'Pin_TrI0gain', 'Pin_TrPD','Pin_TrPDgain',
                            'PresetTime', 'monoE', 'pin_ccd_center_x_pixel','pin_ccd_center_y_pixel',
                            'pin_ccd_tilt_x', 'pin_ccd_tilt_y', 'wavelength', 'waxs_ccd_center_x', 'waxs_ccd_center_y',
                            'waxs_ccd_tilt_x', 'waxs_ccd_tilt_y', 'waxs_ccd_center_x_pixel', 'waxs_ccd_center_y_pixel',
                            'scaler_freq', 'StartTime',                     
                        ]        
            metadata_group = f['/entry/Metadata']
            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            sample_group = f['entry/sample']
            sample_dict = read_group_to_dict(sample_group)
            control_group = f['/entry/control']
            control_dict = read_group_to_dict(control_group)
            Sample["RawData"]["instrument"] = instrument_dict
            Sample["RawData"]["metadata"] = metadata_dict
            Sample["RawData"]["sample"].update(sample_dict)
            Sample["RawData"]["control"] = control_dict
 
        return Sample

def _get_h5_value(h5file, path):
    """Return dataset value at path or None if missing."""
    if path in h5file:
        return h5file[path][()]
    return None

def save_dict_to_hdf5(dic, location, h5file):
    """
    Save a dictionary to an HDF5 file.

    Parameters:
    dic (dict): The dictionary to save.
    filename (str): The name of the HDF5 file.
    """
    def recursively_save_dict_contents_to_group(h5file, path, dic):
        for key, item in dic.items():
            if isinstance(item, dict):
                # Create a new group for nested dictionaries
                logging.debug(f"Creating group: {path} + {key}")
                h5file.create_group(path + key)
                recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
            else:
                # Save numpy arrays and other data types
                h5file[path + key] = item

    recursively_save_dict_contents_to_group(h5file, location, dic)

# # Example usage
# data_dict = {
#     'array': np.array([1, 2, 3]),
#     'value': 42,
#     'nested': {
#         'string': 'hello',
#         'array2': np.array([4, 5, 6])
#     }
# }

#save_dict_to_hdf5(data_dict, 'data.h5')

def load_dict_from_hdf5(hdf_file, location):
    """
    Load a dictionary from an HDF5 file.

    Parameters:
    filename (str): The name of the HDF5 file.

    Returns:
    dict: The loaded dictionary.

    location (str): The location in the HDF5 file to load from.
    'root:DisplayData/'
    """
    def recursively_load_dict_contents_from_group(h5file, path):
        ans = {}
        for key, item in h5file[path].items():
            if isinstance(item, h5py._hl.group.Group):
                ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
            else:
                tempItem = item[()]
                if isinstance(tempItem, bytes):             # Convert bytes to string
                    tempItem = tempItem.decode('utf-8')
                ans[key] = tempItem
        return ans

    return recursively_load_dict_contents_from_group(hdf_file,location)

# Function to recursively read a group and store its datasets in a dictionary
def read_group_to_dict(group):
    data_dict = {}
    for key, item in group.items():
        if isinstance(item, h5py.Dataset):
            # Read the dataset
            data = item[()]
             # Check if the dataset is bytes
            if isinstance(data, bytes):
                # Decode bytes to string
                data = data.decode('utf-8')
            # Check if the dataset is an array with a single element
            elif hasattr(data, 'size') and data.size == 1:
                # Convert to a scalar (number or string)
                data = data.item()
                if isinstance(data, bytes):
                    # Decode bytes to string, the above does not seem to catch this? 
                    data = data.decode('utf-8')
            data_dict[key] = data
        elif isinstance(item, h5py.Group):
            # If the item is a group, recursively read its contents
            data_dict[key] = read_group_to_dict(item)
    return data_dict


# this should not fail if keys on the list are not present
def filter_nested_dict(d, keys_to_keep):
    if isinstance(d, dict):
        return {k: filter_nested_dict(v, keys_to_keep) for k, v in d.items() if k in keys_to_keep and k in d}
    elif isinstance(d, list):
        return [filter_nested_dict(item, keys_to_keep) for item in d]
    else:
        return d    


# def find_NXcanSAS_entries(group, path=''):
#     nxcanSAS_entries = []
    
#     for name, item in group.items():
#         current_path = f"{path}/{name}" if path else name
        
#         # Check if the item is a group
#         if isinstance(item, h5py.Group):
#             # Check if the group has the attribute "NXcanSAS"
#             if 'canSAS_class' in item.attrs:
#                 if(item.attrs['canSAS_class'] == 'SASentry'):
#                     if "definition" in item:
#                         definition_data = item["definition"][()]
#                         # Check if "NXcanSAS" is in the definition data
#                         if isinstance(definition_data, bytes):
#                             definition_data = definition_data.decode('utf-8')
                        
#                         print(f"Definition data: {definition_data}")
#                         if definition_data == 'NXcanSAS':
#                             print(f"Found NXcanSAS entry at: {current_path}")
#                             nxcanSAS_entries.append(current_path)
            
#             # Recursively search within the group
#             nxcanSAS_entries.extend(find_NXcanSAS_entries(item, current_path))
    
#     return nxcanSAS_entries

# this code can find any group which contains listed attributes:values and items:values (strings and variables)
# this is general purpose code for HDF5 - Nexus evaluation
def find_matching_groups(hdf5_file, required_attributes, required_items):
    def check_group(name, obj):
        if isinstance(obj, h5py.Group):
            # Check attributes
            attributes_match = all(
                attr in obj.attrs and obj.attrs[attr] == value
                for attr, value in required_attributes.items()
            )
            
            # Check items
            items_match = True
            for item, expected_value in required_items.items():
                if item in obj:
                    actual_value = obj[item][()]
                    # Unwrap shape-(1,) arrays (Igor/Irena writes definition as 1-element array)
                    if hasattr(actual_value, 'ndim') and actual_value.ndim >= 1 and actual_value.size == 1:
                        actual_value = actual_value.flat[0]
                    if isinstance(actual_value, bytes):
                        actual_value = actual_value.decode('utf-8')
                    if actual_value != expected_value:
                        items_match = False
                        break
                else:
                    items_match = False
                    break
            
            if attributes_match and items_match:
                matching_group_paths.append(name)

    matching_group_paths = []

    hdf5_file.visititems(check_group)

    return matching_group_paths
