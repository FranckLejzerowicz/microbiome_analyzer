# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import dirname, isdir, isfile, islink


def io_update(
        self,
        i_f=None,
        i_d=None,
        o_f=None,
        o_d=None,
        key=None
):
    if not isinstance(key, tuple):
        key = (key,)
    for (IO_fd, val) in [
        (('I', 'f'), i_f),
        (('I', 'd'), i_d),
        (('O', 'f'), o_f),
        (('O', 'd'), o_d)
    ]:
        if key not in self.ios:
            self.ios[key] = {}
        # if IO_fd not in self.ios:
        #     self.ios[IO_fd] = {}
        if not val:
            continue
        if isinstance(val, list):
            self.ios[key].setdefault(IO_fd, set()).update(val)
        elif isinstance(val, str):
            self.ios[key].setdefault(IO_fd, set()).add(val)
        # if isinstance(val, list):
        #     self.ios[IO_fd].setdefault(key, set()).update(val)
        # elif isinstance(val, str):
        #     self.ios[IO_fd].setdefault(key, set()).add(val)


def to_do(
        file: str = None,
        folder: str = None
) -> bool:
    if file:
        file = file.replace('${SCRATCH_FOLDER}', '')
        if isfile(file) or islink(file):
            return False
    if folder:
        folder = folder.replace('${SCRATCH_FOLDER}', '')
        if isdir(folder) or islink(folder):
            return False
    return True


def rep(string: str) -> str:
    return string.replace('${SCRATCH_FOLDER}', '')


def get_roundtrip(io) -> dict:
    roundtrip = {'to': inputs_to_scratch(io), 'from': outputs_back(io)}
    return roundtrip


def inputs_to_scratch(io) -> list:
    rsyncs, mkdirs = set(), set()
    # folders
    if ('I', 'd') in io:
        for folder_ in io[('I', 'd')]:
            folder = folder_.rstrip('/')
            src = folder_.rstrip('/').replace('${SCRATCH_FOLDER}', '')
            mkdirs.add('mkdir -p %s' % folder)
            rsyncs.add('rsync -aqruv %s/ %s' % (src, folder))
    # folders
    if ('O', 'd') in io:
        for folder in io[('O', 'd')]:
            mkdirs.add('mkdir -p %s' % folder.rstrip('/'))
    if ('O', 'f') in io:
        for fil in io[('O', 'f')]:
            mkdirs.add('mkdir -p %s' % dirname(fil))
    # files
    if ('I', 'f') in io:
        for file in io[('I', 'f')]:
            folder = dirname(file)
            src = file.replace('${SCRATCH_FOLDER}', '')
            mkdirs.add('mkdir -p %s' % folder)
            rsyncs.add('rsync -aqruv %s %s' % (src, file))
    return sorted(mkdirs) + sorted(rsyncs)


def outputs_back(io) -> list:
    outbound = set()
    if ('O', 'd') in io:
        # folders
        for folder_ in io[('O', 'd')]:
            folder = folder_.rstrip('/')
            src = folder_.rstrip('/').replace('${SCRATCH_FOLDER}', '')
            cmd = 'mkdir -p %s; rsync -aqruv %s/ %s' % (src, folder, src)
            cmd = 'if [ -d %s ]; then %s; fi' % (folder, cmd)
            outbound.add(cmd)
    if ('O', 'f') in io:
        # files
        for file in io[('O', 'f')]:
            src = file.replace('${SCRATCH_FOLDER}', '')
            folder = dirname(src)
            cmd = 'mkdir -p %s; rsync -aqruv %s %s' % (folder, file, src)
            cmd = 'if [ -f %s ]; then %s; fi' % (file, cmd)
            outbound.add(cmd)
    return sorted(outbound)
