"""Structured error helpers for the pyirena control API.

Every public function in pyirena.api.control returns a plain dict.  When
something goes wrong the dict has an "error" key rather than raising an
exception, so an LLM agent can read the error and decide what to do next.
"""
from __future__ import annotations


def make_error(message: str, suggestion: str = "", code: str = "ERROR") -> dict:
    return {"error": message, "suggestion": suggestion, "code": code}


def no_session(session_id: str) -> dict:
    return make_error(
        f"No session found with id '{session_id}'.",
        suggestion="Call open_dataset() to create a session first.",
        code="NO_SESSION",
    )


def no_model(session_id: str) -> dict:
    return make_error(
        f"Session '{session_id}' has no model selected.",
        suggestion="Call select_model() to choose a fitting model.",
        code="NO_MODEL",
    )


def no_fit(session_id: str) -> dict:
    return make_error(
        f"Session '{session_id}' has no fit results yet.",
        suggestion="Call run_fit() first.",
        code="NO_FIT",
    )


def bad_param(param_name: str, model_name: str) -> dict:
    return make_error(
        f"Parameter '{param_name}' not found in model '{model_name}'.",
        suggestion=(
            "Call get_model_parameters() to see the list of valid parameter names."
        ),
        code="BAD_PARAM",
    )
