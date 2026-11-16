import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('terra_check', ROOT / 'tools/check_terra_wdl.py')
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)


class TerraScopeTests(unittest.TestCase):
    def check_source(self, source):
        import WDL
        return checker.workflow_writers(WDL.parse_document(source))

    def test_workflow_declaration_rejected(self):
        self.assertEqual(len(self.check_source('version 1.0\nworkflow w { File f = write_lines(["x"]) }')), 1)

    def test_workflow_scatter_and_call_input_rejected(self):
        source = '''version 1.0
        task t { input { File f } command { cat ~{f} } }
        workflow w { scatter (i in range(2)) { call t { input: f=write_lines(["x"]) } } }
        '''
        self.assertEqual(len(self.check_source(source)), 1)

    def test_nested_conditional_output_rejected(self):
        source = '''version 1.0
        workflow w { if (true) { File f = write_json({"x": ["y"]}) }
          output { File g = write_lines(["z"]) } }
        '''
        self.assertEqual(len(self.check_source(source)), 2)

    def test_task_writer_allowed(self):
        source = '''version 1.0
        task t { File f = write_lines(["x"]) command { cat ~{f} } }
        workflow w { call t }
        '''
        self.assertEqual(self.check_source(source), [])

    def test_comments_and_literals_do_not_trigger(self):
        source = '''version 1.0
        workflow w { # write_lines(["x"])
          String text = "write_json(x)" }
        '''
        self.assertEqual(self.check_source(source), [])


if __name__ == '__main__':
    unittest.main()
