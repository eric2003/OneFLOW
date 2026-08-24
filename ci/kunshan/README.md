# Kunshan regression workflow

The workflow in `.github/workflows/kunshan-regression.yml` is independent of
the existing Linux, Windows, and documentation workflows. It is intentionally
manual-only because it consumes Slurm resources and repository secrets are not
available to workflows from untrusted forks. It also stops before SSH setup
unless the selected workflow ref is the trusted `master`, `main`, or `hexin`
branch. This repository currently uses `master` as its default branch.

## Required GitHub Actions secrets

| Secret | Purpose |
|---|---|
| `KUNSHAN_SSH_HOST` | SSH host name |
| `KUNSHAN_SSH_PORT` | SSH port |
| `KUNSHAN_SSH_USER` | SSH user |
| `KUNSHAN_SSH_PRIVATE_KEY` | Private key accepted by the cluster |
| `KUNSHAN_SSH_KNOWN_HOSTS` | Pinned OpenSSH known-hosts line |
| `KUNSHAN_REMOTE_ROOT` | Dedicated remote CI root; runs are placed below `runs/` |
| `KUNSHAN_CI_CONFIG` | Absolute path to the private remote configuration file |

If any secret is absent, the private key is invalid, BatchMode SSH login
fails, Slurm commands are unavailable, or the remote configuration is
unreadable, the workflow stops before creating a run directory or submitting
a job and prints a specific Actions error.

## Remote configuration

Copy `config.example.sh` to a private path on the Kunshan shared filesystem,
fill in the target-cluster modules and dependency locations, and set
`KUNSHAN_CI_CONFIG` to that path.

The MPI fixture is kept outside Git because its partitioned grid is large.
`KUNSHAN_MPI4_CASE_ROOT` must point to a complete case directory containing:

- the four-zone grid;
- scripts configured to use that grid;
- an `autotest/` directory containing the accepted four-rank CPU baseline.

Serial and MPI baselines are deliberately separate.

## Execution

Run **OneFLOW Kunshan Regression** from the Actions page and select:

- `cpu-serial`: the three cases in `test/suites/cpu-serial.txt`;
- `mpi4`: the private four-zone/four-rank M6 fixture;
- `all`: serial first, then MPI.

Each workflow run uses:

```text
KUNSHAN_REMOTE_ROOT/runs/<github-run-id>_<attempt>/
```

The source checkout, build, work directories, Slurm logs, accounting output,
and summary remain isolated under that directory. Logs and summaries are also
uploaded as a GitHub Actions artifact for 14 days.
