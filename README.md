# Summer Project 2026

This repo is connected to a Google Drive data folder through a local manifest:

- Drive folder: `Summer 2026_ AD_Research`
- Local dataset path expected by code: `AD analysis/AD datasets/`
- Optional override: set `AD_DATA_DIR` to a synced `AD datasets` folder.

The large data folders are intentionally ignored by git. To check whether the Drive data is available locally, run:

```bash
python "AD analysis/drive_connection/check_data_connection.py"
```

If Google Drive for desktop is installed, make the Drive folder available locally, then either keep/sync it at `AD analysis/AD datasets/` or point the project at the synced folder:

```bash
export AD_DATA_DIR="/path/to/Summer 2026_ AD_Research/AD datasets"
```

You can also create the ignored local links automatically after Drive for desktop exposes the folder:

```bash
python "AD analysis/drive_connection/connect_drive_data.py"
```
