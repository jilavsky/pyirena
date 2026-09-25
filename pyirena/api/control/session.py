"""In-memory session registry for the pyirena control API.

A Session holds one loaded dataset and its associated fitting state.
Sessions live for the lifetime of the calling process — there is no
persistence across restarts.
"""
from __future__ import annotations

import uuid
from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np


@dataclass
class Session:
    session_id: str
    # None for a session created from arrays (open_dataset_from_data): there
    # is no file behind it, so save_* needs an explicit output_path and the
    # summaries report the file as null rather than inventing a name.
    file_path: Optional[str]
    q: np.ndarray
    intensity: np.ndarray
    error: Optional[np.ndarray]
    label: str = ""

    # Q resolution, when the caller supplied it. Stored for provenance and
    # for future resolution smearing; no fitting tool reads it yet.
    dq: Optional[np.ndarray] = None

    # Model state
    model_name: Optional[str] = None
    model: Optional[Any] = None   # UnifiedFitModel | ... depending on model_name

    # Q-range restriction (None = use full data range)
    fit_q_min: Optional[float] = None
    fit_q_max: Optional[float] = None

    # Slit smearing: is_slit_smeared reflects the loaded data; slit_length (1/Å)
    # is applied to the model when a model is created for this session.
    is_slit_smeared: bool = False
    slit_length: float = 0.0

    # Last fit output. A plain dict for most tools, but Modeling and the
    # Carbon model store their own result dataclass here, so anything
    # reading it must handle both.
    last_fit_result: Optional[Any] = None


_SESSIONS: Dict[str, Session] = {}


def create_session(
    file_path: Optional[str],
    q: np.ndarray,
    intensity: np.ndarray,
    error: Optional[np.ndarray] = None,
    label: str = "",
    is_slit_smeared: bool = False,
    slit_length: float = 0.0,
    dq: Optional[np.ndarray] = None,
) -> Session:
    sid = str(uuid.uuid4())[:8]
    session = Session(
        session_id=sid,
        file_path=file_path,
        q=np.asarray(q, dtype=float),
        intensity=np.asarray(intensity, dtype=float),
        error=np.asarray(error, dtype=float) if error is not None else None,
        label=label,
        dq=np.asarray(dq, dtype=float) if dq is not None else None,
        is_slit_smeared=is_slit_smeared,
        slit_length=slit_length,
    )
    _SESSIONS[sid] = session
    return session


def get_session(session_id: str) -> Optional[Session]:
    return _SESSIONS.get(session_id)


def drop_session(session_id: str) -> bool:
    return _SESSIONS.pop(session_id, None) is not None


def all_sessions() -> list[Session]:
    return list(_SESSIONS.values())


def fit_mask(session: Session) -> np.ndarray:
    """Boolean mask selecting the Q points used for fitting."""
    mask = np.ones(len(session.q), dtype=bool)
    if session.fit_q_min is not None:
        mask &= session.q >= session.fit_q_min
    if session.fit_q_max is not None:
        mask &= session.q <= session.fit_q_max
    return mask
