"""IsoTrack calibration handler: one HTCondor job per era."""
import os

from base import HandlerBase


class IsoTrackCalibHandler(HandlerBase):
    task = "isotrack-calib-production"
    script = "automator/scripts/run-isotrack-calib.sh"

    def __init__(self, era, tag, file_list=None, **kwargs):
        super().__init__(**kwargs)
        self.era = era
        self.tag = tag
        self.file_list = file_list

    def inputFiles(self):
        return [self.file_list]

    def getWorkItems(self):
        if not self.file_list or not os.path.isfile(self.file_list):
            raise RuntimeError("file list not found: %s" % self.file_list)
        # The file list is transferred with the job, so only its name is passed.
        args = "%s %s %s" % (self.era, os.path.basename(self.file_list), self.tag)
        return [{"arguments": args}], None
