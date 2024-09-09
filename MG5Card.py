import os
import glob
import re
import tempfile
import subprocess
import shutil

class MG5Card:

    _sub_pattern = re.compile(r'\$SUB:([a-zA-Z0-9_]*)')
    _width_pattern = re.compile(r'=== Results Summary for [^\n]* ===\s+'
                               r'Width\s*:\s*([0-9.Ee+-]+)\s*\+\-\s*([0-9.Ee+-]+)\s*(\S+)\s+'
                               r'Nb of events\s*:\s*(\d+)')

    def __init__(self, path):
        self.path = path
        with open(self.path) as file:
            self.content = file.read()
        self.vars = self._find_vars()

    def _find_vars(self):
        return sorted(set(self._sub_pattern.findall(self.content)))

    def sub(self, vals):
        content = self.content
        for var in self.vars:
            content = content.replace(f'$SUB:{var}', str(vals[var]))
        return content

    def run(self, vals, tmp_workdir=True, capture_output=True):
        if tmp_workdir: vals['workdir'] = tempfile.TemporaryDirectory().name
        content = self.sub(vals)
        datpath = tempfile.mktemp()
        with open(datpath, 'w') as datfile: datfile.write(content)
        proc = subprocess.run(['mg5_aMC', datpath], capture_output=capture_output)
        if tmp_workdir: shutil.rmtree(vals['workdir'])
        os.remove(datpath)
        return proc

    def run_width(self, vals):
        proc = self.run(vals)
        output = proc.stdout.decode()
        width, width_unc, width_unit, nevent = self._width_pattern.search(output).groups()
        return float(width), float(width_unc), width_unit, int(nevent)

def load_cards(basedir='cards'):
    cards = { }
    for path in glob.glob(os.path.join(basedir, '*.dat')):
        cards[os.path.basename(path)] = MG5Card(path)
    return cards

if __name__ == '__main__':
    for path, card in load_cards().items():
        print(path, card.vars, card.content, sep='\n' + '-' * 40 + '\n')
