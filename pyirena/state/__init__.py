"""
State management for pyIrena.

This module handles saving and loading application state across sessions.
"""

from .state_manager import StateManager, get_default_state_file

__all__ = ['StateManager', 'get_default_state_file']
