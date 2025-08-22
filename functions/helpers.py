import sys, fileinput
from pathlib import Path
import datetime as dt
import numpy as np
import pandas as pd
from oslo_concurrency import lockutils
import logging

lockutils.set_defaults(lock_path=Path.home() / "temp")
log = logging.getLogger(__name__)


@lockutils.synchronized('update_l2b_queue', external=True, fair=True)
def update_l2b_queue(image, flag, file_queue):
    image = str(image)
    for line in fileinput.input(file_queue, inplace=True):
        if image in line:
            line = f"{image} {flag}\n"

        sys.stdout.write(line)


def read_l2b_filequeue(l2b_filequeue):
    l2a_files = []
    status = []
    with open(l2b_filequeue, "r") as fid:
        for line in fid:
            l2a_file, mode = line.rstrip("\n\r").split(" ")
            l2a_files.append(l2a_file)
            status.append(mode)

    return l2a_files, status


def write_l2b_filequeue(l2b_filequeue, provenances_to_process, region="europe"):

    l2a_files = []
    modes = []
    for provenance_file in provenances_to_process:
        temp_list = pd.read_csv(provenance_file)
        temp_list = temp_list[["file", "mode"]]
        for i, (file, mode) in temp_list.iterrows():
            if "BOA" in file and region in file:
                l2a_files.append(file)
                modes.append(mode)

    l2a_files_set = set(l2a_files)
    l2a_files = np.array(l2a_files)
    modes = np.array(modes)
    existing_l2a_files, _ = read_l2b_filequeue(l2b_filequeue)
    for l2a_file in l2a_files_set:
        valid = l2a_file == l2a_files
        np.sum(valid)
        if "create" in modes[valid]:
            flag = "NEW"
        else:
            flag = "UPDATED"

        if l2a_file in existing_l2a_files:
            update_l2b_queue(l2a_file, flag, l2b_filequeue)
        else:
            with open(l2b_filequeue, "a") as fid:
                fid.write(f"{l2a_file} {flag}\n")
                fid.flush()


def create_l2b_filequeue(l2a_storage, l2b_storage, l2b_filequeue):
    l2a_storage = Path(l2a_storage)
    l2b_storage = Path(l2b_storage)
    l2b_filequeue = Path(l2b_filequeue)

    with open(l2b_filequeue, "w") as fid:
        tiles = l2b_storage.glob("X*")
        for tile in tiles:
            l2b_files = tile.glob("*_LAI.tif")
            for l2b_file in l2b_files:
                datestr, level, satellite, _ = l2b_file.stem.split("_")
                l2a_file = l2a_storage / tile.stem / \
                           f"{datestr}_{level}_{satellite}_BOA.tif"

                if l2a_file.exists():
                    fid.write(f"{l2a_file} DONE\n")

        fid.flush()



def get_provenances_to_process(provenance_folder,
                               processed_provenance_filelist):

    provenance_folder = Path(provenance_folder)
    provenance_files = provenance_folder.glob("*.csv")
    processed_provenance = []
    with open(processed_provenance_filelist, "r") as fid:
        for line in fid:
            processed_provenance.append(line.rstrip("\n\r"))

    to_process = [file for file in provenance_files
                  if str(file) not in processed_provenance]

    return to_process


def update_processed_provenance(processed_provenance_filelist,
                                provenances_to_process):
    with open(processed_provenance_filelist, "a") as fid:
        for provenance in provenances_to_process:
            fid.write(f"{provenance}\n")

        fid.flush()


def _check_bracketting_date(date_str, start_date=None, end_date=None):
    skip = False
    date_obj = dt.datetime.strptime(date_str[:8], "%Y%m%d").date()
    if not isinstance(start_date, type(None)):
        if date_obj < start_date:
            log.info(f"Acquisition date {date_obj} is before {start_date}")
            skip = True
    if not isinstance(end_date, type(None)):
        if date_obj > end_date:
            log.info(f"Acquisition date {date_obj} is after {end_date}")
            skip = True

    return skip
