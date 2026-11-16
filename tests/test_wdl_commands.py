"""Render every task after cloud-to-local File mapping and inspect its CLI."""
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest

import WDL

ROOT = pathlib.Path(__file__).resolve().parents[1]


class LocalTaskStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename


class LocalizedCommandTests(unittest.TestCase):
    def test_all_task_commands_localize_and_quote_files(self):
        document = WDL.load(str(ROOT / 'workflows/xyalign.wdl'))
        for task in document.tasks:
            for optional_present in (False, True):
                with self.subTest(task=task.name, optional_present=optional_present), tempfile.TemporaryDirectory() as directory:
                    root = pathlib.Path(directory)
                    values = dict(cpus=2, memory_gb=16, disk_gb=100, preemptible_attempts=0,
                                  docker_image='image:sha', sample_id='sample',
                                  x_chromosome='chrX', y_chromosome='chrY', y_status='present',
                                  remapping_requested=False, characterize=False, call_variants=False,
                                  single_end=False, java_memory='8g', window_size=100,
                                  num_permutations=10, num_bootstraps=10, chromosomes=['chr19', 'chrX', 'chrY'])
                    mapping = {}
                    for binding in task.available_inputs:
                        declaration = binding.value
                        if isinstance(declaration.type, WDL.Type.File):
                            if declaration.type.optional and not optional_present:
                                continue
                            cloud = 'gs://input-bucket/' + binding.name
                            local = root / (binding.name + " ' $(touch PWNED)")
                            local.write_text('localized input\n')
                            values[binding.name] = cloud
                            mapping[cloud] = str(local)
                    values = {k: v for k, v in values.items() if k in task.available_inputs}
                    environment = WDL.values_from_json(values, task.available_inputs)
                    environment = WDL.Value.rewrite_env_paths(environment, lambda value: mapping[value.value])
                    # Set absent optional values exactly as a task runtime does.
                    for declaration in task.inputs:
                        if declaration.name not in environment:
                            environment = environment.bind(declaration.name, WDL.Value.Null())
                    stdlib = LocalTaskStdLib('1.0', write_dir=str(root))
                    for declaration in task.postinputs:
                        environment = environment.bind(declaration.name, declaration.expr.eval(environment, stdlib))
                    command = task.command.eval(environment, stdlib).value
                    capture = root / 'capture.py'
                    capture.write_text('import json,sys\nopen("arguments.json", "w").write(json.dumps(sys.argv[1:]))\n')
                    command = command.replace('python /opt/sexchromosome/scripts/xyalign_runner.py',
                                              sys.executable + ' ' + str(capture))
                    result = subprocess.run(['bash', '-c', command], cwd=directory, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    arguments = json.loads((root / 'arguments.json').read_text())
                    self.assertNotIn('gs://', ' '.join(arguments))
                    self.assertFalse((root / 'PWNED').exists())
                    for cloud, local in mapping.items():
                        if cloud.endswith('/analysis_log'):
                            continue
                        self.assertIn(local, arguments)
                    self.assertEqual('--reference-index' in arguments, optional_present)
                    if task.name == 'PrepareXYReferences':
                        self.assertEqual('--reference-mask' in arguments, optional_present)

    def test_struct_files_retain_file_types(self):
        document = WDL.load(str(ROOT / 'workflows/xyalign.wdl'))
        values = {'XYalign.samples': [{'sample_id': 's', 'bam': 'gs://bucket/s.bam',
                                     'bam_index': 'gs://bucket/s.bam.bai'}],
                  'XYalign.reference_fasta': 'gs://bucket/ref.fa'}
        environment = WDL.values_from_json(values, document.workflow.available_inputs, namespace='XYalign')
        sample = environment['samples'].value[0]
        self.assertIsInstance(sample.value['bam'], WDL.Value.File)
        self.assertIsInstance(sample.value['bam_index'], WDL.Value.File)
        localized = WDL.Value.rewrite_env_paths(environment, lambda value: '/localized/' + value.value.rsplit('/', 1)[-1])
        self.assertEqual(localized['samples'].value[0].value['bam'].value, '/localized/s.bam')


if __name__ == '__main__':
    unittest.main()
