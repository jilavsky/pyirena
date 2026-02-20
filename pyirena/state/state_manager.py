"""
State manager for pyIrena application.

Handles saving and loading application state in a hierarchical JSON format.
The state file is human-readable and can be edited manually if needed.
"""

import json
import os
from pathlib import Path
from typing import Dict, Any, Optional
from copy import deepcopy


def get_default_state_file() -> Path:
    """
    Get the default state file path.

    Returns the path to ~/.pyirena/state.json
    """
    home = Path.home()
    state_dir = home / '.pyirena'
    state_dir.mkdir(exist_ok=True)
    return state_dir / 'state.json'


class StateManager:
    """
    Manages application state for pyIrena.

    The state is stored in a hierarchical JSON structure:
    {
        "version": "1.0",
        "unified_fit": { ... },
        "other_tool": { ... }
    }

    Each tool can have its own section in the state file.
    """

    # Default state for the entire application
    DEFAULT_STATE = {
        "version": "1.0",
        "data_selector": {
            "last_folder": "",
            "error_fraction": 0.05,   # uncertainty = I × error_fraction when file has no error column
        },
        "sizes": {
            # schema_version is bumped whenever a default value changes so that
            # old saved states can be migrated automatically on load.
            "schema_version": 3,
            "r_min": 10.0,
            "r_max": 1000.0,
            "n_bins": 200,           # was 50 in schema_version 1
            "log_spacing": True,     # was False in schema_version 1
            "shape": "sphere",
            "contrast": 1.0,
            "aspect_ratio": 1.0,           # used when shape == 'spheroid'
            "background": 0.0,
            "error_scale": 1.0,            # new in schema_version 2
            "method": "regularization",    # 'maxent' | 'regularization' | 'tnnls' | 'mcsas'
            "maxent_sky_background": 1e-6,
            "maxent_stability": 0.01,
            "maxent_max_iter": 1000,
            "regularization_evalue": 1.0,
            "regularization_min_ratio": 1e-4,
            "tnnls_approach_param": 0.95,
            "tnnls_max_iter": 1000,
            # McSAS Monte Carlo parameters (new in schema_version 3)
            "mcsas_n_repetitions": 10,
            "mcsas_convergence": 1.0,
            "mcsas_max_iter": 100000,
            # Power-law background: B·q^(-P) + flat background
            "power_law_B": 0.0,
            "power_law_P": 4.0,
            "power_law_q_min": None,   # Q range for fitting B and/or P
            "power_law_q_max": None,
            "background_q_min": None,  # Q range for fitting flat background
            "background_q_max": None,
        },
        "unified_fit": {
            "num_levels": 1,
            "levels": [
                {
                    "level": 1,
                    "G": {
                        "value": 100.0,
                        "fit": False,
                        "low_limit": None,
                        "high_limit": None
                    },
                    "Rg": {
                        "value": 100.0,
                        "fit": False,
                        "low_limit": None,
                        "high_limit": None
                    },
                    "B": {
                        "value": 0.01,
                        "fit": False,
                        "low_limit": None,
                        "high_limit": None
                    },
                    "P": {
                        "value": 4.0,
                        "fit": False,
                        "low_limit": None,
                        "high_limit": None
                    },
                    "ETA": {
                        "value": 0.0,
                        "fit": False,
                        "low_limit": None,
                        "high_limit": None
                    },
                    "PACK": {
                        "value": 0.0,
                        "fit": False,
                        "low_limit": None,
                        "high_limit": None
                    },
                    "RgCutoff": 0.0,
                    "correlated": False,
                    "estimate_B": False
                },
                # Levels 2-5 with same structure
                {
                    "level": 2,
                    "G": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "Rg": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "B": {"value": 0.01, "fit": False, "low_limit": None, "high_limit": None},
                    "P": {"value": 4.0, "fit": False, "low_limit": None, "high_limit": None},
                    "ETA": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "PACK": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "RgCutoff": 0.0,
                    "correlated": False,
                    "estimate_B": False
                },
                {
                    "level": 3,
                    "G": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "Rg": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "B": {"value": 0.01, "fit": False, "low_limit": None, "high_limit": None},
                    "P": {"value": 4.0, "fit": False, "low_limit": None, "high_limit": None},
                    "ETA": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "PACK": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "RgCutoff": 0.0,
                    "correlated": False,
                    "estimate_B": False
                },
                {
                    "level": 4,
                    "G": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "Rg": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "B": {"value": 0.01, "fit": False, "low_limit": None, "high_limit": None},
                    "P": {"value": 4.0, "fit": False, "low_limit": None, "high_limit": None},
                    "ETA": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "PACK": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "RgCutoff": 0.0,
                    "correlated": False,
                    "estimate_B": False
                },
                {
                    "level": 5,
                    "G": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "Rg": {"value": 100.0, "fit": False, "low_limit": None, "high_limit": None},
                    "B": {"value": 0.01, "fit": False, "low_limit": None, "high_limit": None},
                    "P": {"value": 4.0, "fit": False, "low_limit": None, "high_limit": None},
                    "ETA": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "PACK": {"value": 0.0, "fit": False, "low_limit": None, "high_limit": None},
                    "RgCutoff": 0.0,
                    "correlated": False,
                    "estimate_B": False
                }
            ],
            "background": {
                "value": 1e-6,
                "fit": False
            },
            "cursor_left": None,  # Will be set when data is loaded
            "cursor_right": None,  # Will be set when data is loaded
            "update_auto": False,
            "display_local": False,
            "no_limits": False,
            "skip_fit_check": False,
            "store_local": False
        }
    }

    def __init__(self, state_file: Optional[Path] = None):
        """
        Initialize the state manager.

        Args:
            state_file: Path to state file. If None, uses default location.
        """
        self.state_file = state_file or get_default_state_file()
        self.state = deepcopy(self.DEFAULT_STATE)
        self.load()

    def load(self) -> bool:
        """
        Load state from file.

        Returns:
            True if state was loaded, False if using defaults
        """
        if not self.state_file.exists():
            print(f"State file not found: {self.state_file}")
            print("Using default state")
            return False

        try:
            with open(self.state_file, 'r') as f:
                loaded_state = json.load(f)

            # Capture schema versions BEFORE merging so we can detect old files.
            # _merge_state gives loaded values priority over defaults, so after
            # merging the 'schema_version' key would reflect DEFAULT_STATE (not
            # the on-disk value) whenever the key is absent in the loaded file.
            loaded_sizes_version = loaded_state.get('sizes', {}).get('schema_version', 1)

            # Merge loaded state with defaults (in case new fields were added)
            self.state = self._merge_state(self.DEFAULT_STATE, loaded_state)
            # Apply any default-value migrations for changed schema versions
            self._migrate_state(loaded_sizes_version)
            print(f"Loaded state from: {self.state_file}")
            return True

        except Exception as e:
            print(f"Error loading state file: {e}")
            print("Using default state")
            return False

    def save(self) -> bool:
        """
        Save current state to file.

        Returns:
            True if successful, False otherwise
        """
        try:
            # Ensure directory exists
            self.state_file.parent.mkdir(parents=True, exist_ok=True)

            # Write with pretty formatting
            with open(self.state_file, 'w') as f:
                json.dump(self.state, f, indent=2)

            print(f"Saved state to: {self.state_file}")
            return True

        except Exception as e:
            print(f"Error saving state file: {e}")
            return False

    def get(self, tool: str, key: Optional[str] = None, default: Any = None) -> Any:
        """
        Get state for a tool or specific key.

        Args:
            tool: Tool name (e.g., "unified_fit")
            key: Optional key within tool state
            default: Default value if not found

        Returns:
            State value or default
        """
        tool_state = self.state.get(tool, {})

        if key is None:
            return tool_state

        return tool_state.get(key, default)

    def set(self, tool: str, key: str, value: Any):
        """
        Set state for a tool.

        Args:
            tool: Tool name (e.g., "unified_fit")
            key: Key within tool state
            value: Value to set
        """
        if tool not in self.state:
            self.state[tool] = {}

        self.state[tool][key] = value

    def update(self, tool: str, state_dict: Dict[str, Any]):
        """
        Update multiple state values for a tool.

        Args:
            tool: Tool name (e.g., "unified_fit")
            state_dict: Dictionary of key-value pairs to update
        """
        if tool not in self.state:
            self.state[tool] = {}

        self.state[tool].update(state_dict)

    def reset(self, tool: Optional[str] = None):
        """
        Reset state to defaults.

        Args:
            tool: Tool to reset. If None, resets all tools.
        """
        if tool is None:
            self.state = deepcopy(self.DEFAULT_STATE)
        else:
            if tool in self.DEFAULT_STATE:
                self.state[tool] = deepcopy(self.DEFAULT_STATE[tool])

    def export_tool_state(self, tool: str, export_path: Path) -> bool:
        """
        Export tool state to a separate file.

        This is useful for sharing fit parameters or creating presets.

        Args:
            tool: Tool name (e.g., "unified_fit")
            export_path: Path to export file

        Returns:
            True if successful, False otherwise
        """
        try:
            tool_state = self.state.get(tool, {})

            with open(export_path, 'w') as f:
                json.dump(tool_state, f, indent=2)

            print(f"Exported {tool} state to: {export_path}")
            return True

        except Exception as e:
            print(f"Error exporting state: {e}")
            return False

    def import_tool_state(self, tool: str, import_path: Path) -> bool:
        """
        Import tool state from a separate file.

        Args:
            tool: Tool name (e.g., "unified_fit")
            import_path: Path to import file

        Returns:
            True if successful, False otherwise
        """
        try:
            with open(import_path, 'r') as f:
                tool_state = json.load(f)

            self.state[tool] = tool_state
            print(f"Imported {tool} state from: {import_path}")
            return True

        except Exception as e:
            print(f"Error importing state: {e}")
            return False

    def _migrate_state(self, loaded_sizes_version: int = None):
        """
        Upgrade saved state to current schema.

        When default values change (e.g. n_bins 50→200) we bump
        ``schema_version`` in DEFAULT_STATE and reset the affected fields
        here so that users automatically get the new defaults rather than
        having the old on-disk values silently persist.

        Args:
            loaded_sizes_version: The schema_version read from the on-disk file
                *before* merging with DEFAULT_STATE.  If None, falls back to
                reading schema_version from the (already-merged) state dict,
                which is only correct when called outside of load().
        """
        sizes = self.state.get('sizes', {})
        # Use the pre-merge version when available (passed from load()); otherwise
        # fall back to the merged state value (used when called standalone).
        stored_version = loaded_sizes_version if loaded_sizes_version is not None \
            else sizes.get('schema_version', 1)
        target_version = self.DEFAULT_STATE['sizes']['schema_version']

        if stored_version < 2 <= target_version:
            # schema_version 1 → 2: n_bins default 50→200, log_spacing False→True,
            # error_scale field added.
            sizes['n_bins']      = self.DEFAULT_STATE['sizes']['n_bins']
            sizes['log_spacing'] = self.DEFAULT_STATE['sizes']['log_spacing']
            sizes['error_scale'] = self.DEFAULT_STATE['sizes']['error_scale']
            sizes['schema_version'] = 2
            self.state['sizes'] = sizes

        if stored_version < 3 <= target_version:
            # schema_version 2 → 3: McSAS Monte Carlo parameters added.
            sizes['mcsas_n_repetitions']   = self.DEFAULT_STATE['sizes']['mcsas_n_repetitions']
            sizes['mcsas_convergence']     = self.DEFAULT_STATE['sizes']['mcsas_convergence']
            sizes['mcsas_max_iter']        = self.DEFAULT_STATE['sizes']['mcsas_max_iter']
            sizes['schema_version'] = 3
            self.state['sizes'] = sizes

    def _merge_state(self, default: Dict, loaded: Dict) -> Dict:
        """
        Merge loaded state with default state.

        This ensures that new fields in DEFAULT_STATE are present even
        if they weren't in the loaded state file.

        Args:
            default: Default state dictionary
            loaded: Loaded state dictionary

        Returns:
            Merged state dictionary
        """
        merged = deepcopy(default)

        for key, value in loaded.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = self._merge_state(merged[key], value)
            else:
                merged[key] = value

        return merged
