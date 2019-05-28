#!/bin/bash
#
# Copyright(c) 2019 Intel Corporation
# SPDX - License - Identifier: BSD - 2 - Clause - Patent
#
set -x

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$script_dir"

call_help(){
    cat <<EOF
Usage: ./.cleanup.sh [OPTION] .. -- [FILES|DIRECTORIES]
    -a, --all       Runs cleanup on all files in tree
    -b <branch>, --branch=<branch>  Use a different branch as origin [master]
    -c, --commit    Commits changes as "Cleanup: Remove trailing whitespace"
    -d, --dry-run   Prints out the files that would have been modified.
                        Does not apply if FILES or DIRECTORIES are present
    -f, --force     Cleans up files regardles of commit status
    -h, --help      This text
    -s, --spaces    Squashes multiple blank lines into one line (cat -s)

    All files and directories are relative to the root of the repository
EOF
}

print_mesage(){
    local FILE="stderr"
    [[ "$1" == "stdout" ]] &&
        FILE=stdout &&
            shift
    # Remove leading tabs from the message, if tabs are desired, use cat
    local message="$(sed -e 's/^[ \t]*//' <<< ${1:-Unknown Message})"
    if [[ $FILE == stdout ]]; then
        echo -e "$message"
    else
        echo -e "$message" >&2
    fi
}

die(){
    local message="$(sed -e 's/^[ \t]*//' <<< ${1:-Unknown error})"
    local exit_code="${2:-1}"
    print_mesage "$message"
    exit $exit_code
}

while true; do
    case $1 in
        --help|-h )
            call_help
            exit ;;
        --all|-a )
            clean_all=y
            shift ;;
        --branch=*)
            branch="${1#*=}"
            print_mesage "Using $branch as origin"
            shift ;;
        -b )
            [[ -z "$2" ]] && die "No branch specified"
            branch="$2"
            print_mesage "Using $branch as origin"
            shift 1;;
        --commit|-c )
            commit_cleanup=y
            shift ;;
        --force|-f )
            force=y
            shift ;;
        --dry-run|-d )
            dry_run=y
            shift ;;
        --spaces|-s)
            squash_space=you
            shift
            ;;
        -- )
            shift
            break ;;
        -* )
            print_mesage "Error, unknown option: '$1'."
            exit 1 ;;
        * )
            break ;;
    esac
done

get_git_folder() {
    local repo_folder="$PWD"
    local base_folder="$(basename $PWD)"
    if [[ -d ".git" ]]; then
        echo "$repo_folder"
    elif [[ -d "../.git" ]]; then
        echo "${repo_folder%$base_folder}"
    else
        die "Can't find .git in current folder nor in parent directory"
    fi
}

get_diff_files() {
    local ovc_name=$(git remote -v | grep -i openvisualcloud | head -n1 | cut -f1)
    ovc_name=${ovc_name:-origin}
    local branch=${branch:-master}
    git diff --name-only $branch | grep -iE "*.sh|*.bat"
    for file in $(git diff --name-only $branch | sed '/*.png/d;/*.git/d;/*.vs*/d;/Build*/d;/Bin*/d'); do
        [[ -f "$file" ]] && echo "$file"
    done
}

if ! command -v find &>/dev/null; then
    die "find not found in PATH. This is required for this script"
fi

if [[ $(git status --porcelain -uno | head -c1 | wc -c) -ne 0 ]] &&
    [[ "$force" != "y" ]]; then
    die "Please commit any uncommitted changes before running this script
    If you want to run this command, use --force
    However, if you do use --force, there is no guarentee if sed fails"
fi

if [[ $(find . -name "*.bak" | head -c1 | wc -c) -ne 0 ]] &&
    [[ "$force" != "y" ]]; then
    die "Old sed backup files (*.bak) found. Please either delete them
    or properly clean them up. To have the script clean these up, use --force"
fi

check_directories=()
check_files=()
ignore_names=( -type f \( -name \'*.sh\' -or -name \'*.bat\' -or \
    ! -path \'*.git/*\' ! -path \'*.vs/*\' ! -path \'*.vscode/*\' \
    ! -path \'*Bin/*\' ! -path \'*Build/*\' ! -name \'*.exe\' \
    ! -name \'*.dll\' ! -name \'*.a\' \! -name \'*.so\' ! -name \'*.lib\' \
    ! -name \'*.png\' ! -name \'*.o\' \) )

if [[ "$#" -ne 0 ]]; then
    while [[ $# -ne 1 ]]; do
        for ford in "$@"; do
            if [[ -d "$ford" ]]; then
                check_directories+=( "$ford" )
            elif [[ -f "$ford" ]]; then
                check_files+=( "$ford" )
            else
                print_mesage "$ford was not found"
            fi
        done
        shift
    done
    if [[ -n "$check_directories" ]]; then
        for dir in "${check_directories[@]}"; do
            find "$dir" ${ignore_names[@]} -exec sed -i.bak 's/[[:space:]]*$//' {} +
            [[ "$squash_space" == "y" ]] &&
                find "$dir" ${ignore_names[@]} -exec bash -c 'cat -s <<<"$(<"$0")" >"$0"' '{}' \;
            find . -type f -name "*.bak" -delete
        done
    fi
    if [[ -n "$check_files" ]]; then
        sed -i.bak 's/[[:space:]]*$//' "$check_files"
        for file in "${check_files[@]}"; do
            cat -s <<<"$(<"$file")" >"$file"
        done
        find . -type f -name "*.bak" -delete
    fi
else
    if [[ "$clean_all" == "y" ]]; then
        if [[ "$dry_run" == "y" ]]; then
            find . ${ignore_names[@]} \
                -exec grep -l '[[:blank:]]$' {} +
        else
            find . ${ignore_names[@]} -exec sed -i.bak 's/[[:space:]]*$//' {} + \
                -exec bash -c 'cat -s <<<"$(<"$0")" > "$0" && echo "">> "$0"' '{}' \;
            find . -type f -name "*.bak" -delete
        fi
    elif [[ $(get_diff_files | head -c1 | wc -c) -ne 0 ]]; then
        if [[ "$dry_run" == "y" ]]; then
            grep -l '[[:blank:]]$' $(get_diff_files)
        else
            find . -type f -name "*.bak" -delete
            if ! sed -i.bak 's/[[:space:]]*$//' $(get_diff_files); then
                find . -type f -name "*.bak" -exec bash -c 'mv $0 $(basename "$0" .bak)' '{}' \;
                die "Failed to modify the changed files, please make sure none of them opened
                and you have permission to modify them"
            else
                find . -type f -name "*.bak" -delete
                for file in $(get_diff_files); do
                    cat -s <<<"$(<"$file")" > "$file"
                done
            fi
        fi
    else
        print_mesage "Nothing to do"
    fi
fi

if [[ "$commit_cleanup" == "y" ]]; then
    if [[ $(git status --porcelain -uno | head -c1 | wc -c) -ne 0 ]]; then
        git commit -am "Cleanup: Remove trailing whitespace"
    else
        print_mesage "Nothing to commit"
    fi
fi

# git merge-base --is-ancestor $(git log --follow -1 --format=%H -- "$file") HEAD
# For potential --fixup option
