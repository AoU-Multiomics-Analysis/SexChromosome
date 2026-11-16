import importlib.util
import contextlib
import io
import pathlib
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('runner', ROOT / 'scripts/xyalign_runner.py')
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class LocalizationTests(unittest.TestCase):
    def test_unresolved_cloud_uri_fails(self):
        for value in ['gs://bucket/sample.bam', 's3://bucket/sample.bam', 'https://host/sample.bam']:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'Localization error'):
                runner.require_file(value)

    def test_localized_file_with_shell_characters_uses_fixed_alias(self):
        with tempfile.TemporaryDirectory() as directory:
            source = pathlib.Path(directory) / "sample ' $(touch BAD).bam"
            source.write_bytes(b'bam')
            target = pathlib.Path(directory) / 'input.bam'
            runner.stage_file(str(source), str(target))
            self.assertEqual(target.read_bytes(), b'bam')
            self.assertTrue(target.is_symlink())

    def test_empty_input_rejected(self):
        with tempfile.NamedTemporaryFile() as handle:
            with self.assertRaisesRegex(ValueError, 'empty'):
                runner.require_file(handle.name)

    def test_reject_unsafe_tokens(self):
        for token in ['$(touch BAD)', '../sample', 'a;echo', 'a\nb', 'sample\n', '--option', '.', '..']:
            with self.subTest(token=token), self.assertRaises(ValueError):
                runner.safe_token(token, 'sample ID')

    def test_cloud_to_local_mapping_then_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            local = pathlib.Path(directory) / 'localized' / 'reads.bam'
            local.parent.mkdir()
            local.write_bytes(b'localized BAM')
            mapping = {'gs://input-bucket/reads.bam': str(local)}
            target = pathlib.Path(directory) / 'input.bam'
            runner.stage_file(mapping['gs://input-bucket/reads.bam'], str(target))
            self.assertEqual(target.read_bytes(), b'localized BAM')

    def test_chromosome_list_validated(self):
        with tempfile.NamedTemporaryFile(mode='w+') as handle:
            handle.write('chr19\nchrX\nchrX\nchrY\n')
            handle.flush()
            with self.assertRaisesRegex(ValueError, 'unique'):
                runner.read_chromosomes(handle.name, 'chrX', 'chrY')


class ReadGroupTests(unittest.TestCase):
    def test_reject_ids_that_upstream_cannot_parse(self):
        for value in ['lane:1', 'lane 1', '../lane', 'lane/1', 'lane\\1']:
            with self.subTest(value=value), self.assertRaises(ValueError):
                runner.validate_read_groups([{'ID': value, 'SM': 'sample'}])

    def test_safe_read_group_is_accepted(self):
        runner.validate_read_groups([{'ID': 'lane_1.2-3', 'SM': 'sample', 'PL': 'ILLUMINA'}])


class ReferenceArchiveTests(unittest.TestCase):
    def test_reject_archive_path_escape(self):
        import io
        import tarfile
        with tempfile.TemporaryDirectory() as directory:
            archive = pathlib.Path(directory) / 'bad.tar'
            with tarfile.open(str(archive), 'w') as output:
                member = tarfile.TarInfo('../escape')
                member.size = 1
                output.addfile(member, io.BytesIO(b'x'))
            with self.assertRaisesRegex(ValueError, 'archive'):
                runner.extract_references(str(archive), directory)

    def test_extract_complete_reference_archive(self):
        import tarfile
        with tempfile.TemporaryDirectory() as directory:
            archive = pathlib.Path(directory) / 'refs.tar'
            with tarfile.open(str(archive), 'w') as output:
                for prefix in ('xx', 'xy'):
                    for suffix in [''] + runner.INDEX_SUFFIXES:
                        member = tarfile.TarInfo('references/' + prefix + '.fa' + suffix)
                        member.size = 1
                        output.addfile(member, io.BytesIO(b'x'))
            runner.extract_references(str(archive), directory)
            self.assertEqual((pathlib.Path(directory) / 'references/xx.fa').read_bytes(), b'x')

    def test_remap_requires_explicit_y_status(self):
        parser = runner.make_parser()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parser.parse_args(['remap', '--bam', 'x.bam', '--bam-index', 'x.bai',
                               '--reference', 'ref.fa', '--sample-id', 's',
                               '--prepared-references', 'refs.tar', '--output-dir', 'out'])


if __name__ == '__main__':
    unittest.main()
