#!/bin/sh

DG_SCRIPT_PATH="${BASH_SOURCE[0]:-${(%):-%x}}"
DG_ROOT="$(cd "$(dirname "$DG_SCRIPT_PATH")/.." && pwd -P)"

if [ -n "$DG_OPENFOAM_BASHRC" ]
then
    FOAM_BASHRC="$DG_OPENFOAM_BASHRC"
elif [ -n "$WM_PROJECT_DIR" ] && [ -f "$WM_PROJECT_DIR/etc/bashrc" ]
then
    FOAM_BASHRC="$WM_PROJECT_DIR/etc/bashrc"
elif [ -f "$DG_ROOT/../OpenFOAM-v2412/etc/bashrc" ]
then
    FOAM_BASHRC="$DG_ROOT/../OpenFOAM-v2412/etc/bashrc"
else
    echo "Error: could not locate OpenFOAM bashrc." 1>&2
    echo "Set DG_OPENFOAM_BASHRC or source OpenFOAM manually first." 1>&2
    return 1 2>/dev/null || exit 1
fi

. "$FOAM_BASHRC" WM_COMPILE_OPTION=Debug || return 1 2>/dev/null || exit 1
