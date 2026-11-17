#!/usr/bin/env bash

set -Eeuo pipefail

repo_root="$(git -C "$(dirname "${BASH_SOURCE[0]}")/.." rev-parse --show-toplevel)"
state_root="${XDG_STATE_HOME:-${HOME:?HOME is not set}/.local/state}/oneflow"
mkdir -p "$state_root"

exec 9>"$state_root/master-sync.lock"
if ! flock -n 9; then
    printf '[%s] another master sync is running; skip\n' "$(date -Is)"
    exit 0
fi

log() {
    printf '[%s] %s\n' "$(date -Is)" "$*"
}

die() {
    log "ERROR: $*"
    exit 1
}

dry_run=false
case "${1:-}" in
    "") ;;
    --dry-run) dry_run=true ;;
    --help|-h)
        cat <<'EOF'
Usage: scripts/sync-master.sh [--dry-run]

Fetch origin/master and upstream/master, then fast-forward origin/master when
it is safely behind upstream/master. The script never changes local master or
dev. --dry-run fetches remotes and reports the planned push without pushing.
EOF
        exit 0
        ;;
    *) die "unknown argument: $1" ;;
esac

cd "$repo_root"
origin_url="$(git remote get-url origin)" || die "remote 'origin' is missing"
upstream_url="$(git remote get-url upstream)" || die "remote 'upstream' is missing"
case "$origin_url" in
    *github.com:lql341/OneFLOW.git|*github.com/lql341/OneFLOW.git) ;;
    *) die "unexpected origin URL: $origin_url" ;;
esac
case "$upstream_url" in
    *github.com:eric2003/OneFLOW.git|*github.com/eric2003/OneFLOW.git) ;;
    *) die "unexpected upstream URL: $upstream_url" ;;
esac

log "fetching master from origin and upstream"
git fetch --no-tags origin \
    +refs/heads/master:refs/remotes/origin/master
git fetch --no-tags upstream \
    +refs/heads/master:refs/remotes/upstream/master

origin_master="$(git rev-parse refs/remotes/origin/master)"
upstream_master="$(git rev-parse refs/remotes/upstream/master)"

if ! git merge-base --is-ancestor "$origin_master" "$upstream_master"; then
    die "origin/master is not an ancestor of upstream/master; refusing to overwrite fork-only commits"
fi
log "state: origin=$origin_master upstream=$upstream_master"
if [[ "$origin_master" != "$upstream_master" ]]; then
    if [[ "$dry_run" == true ]]; then
        log "dry-run: would fast-forward origin/master to $upstream_master"
    else
        log "fast-forwarding origin/master"
        git push origin "$upstream_master:refs/heads/master"
    fi
fi

if [[ "$dry_run" == true ]]; then
    log "dry-run complete"
else
    log "remote sync complete; origin/master matches upstream/master"
fi
