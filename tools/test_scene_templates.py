"""Run with python -m unittest discover -s tools -p test_scene_templates.py."""

from pathlib import Path
import tempfile
import unittest

from scene_templates import load_scene


class BodyTemplatesTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='astraxis_templates_')
        self.addCleanup(self.temp.cleanup)
        self.assets = Path(self.temp.name)
        (self.assets / 'scenes').mkdir()
        (self.assets / 'bodies').mkdir()
        self.path = self.assets / 'scenes/test.toml'
        self.template = self.assets / 'bodies/test.toml'
        self.template.write_text('radii_km = [10, 8]\n[atmosphere]\nheight_km = 1\ngain = 2\n', encoding='utf-8')

    def load(self, body):
        self.path.write_text("name = 'Test'\n[[bodies]]\nname = 'A'\n" + body, encoding='utf-8')
        return load_scene(self.path)['bodies']

    def test_overlay_and_instance_isolation(self):
        a, b = self.load("template = 'bodies/test.toml'\nradii_km = [4]\n"
                         "[bodies.atmosphere]\ngain = 0\n"
                         "[[bodies]]\nname = 'B'\ntemplate = 'bodies/test.toml'\n")
        self.assertEqual(a['radii_km'], [4])
        self.assertEqual(a['atmosphere'], {'height_km': 1, 'gain': 0})
        self.assertEqual(b['radii_km'], [10, 8])
        self.assertEqual(b['atmosphere']['gain'], 2)
        self.assertNotIn('template', a)

    def test_unset_then_replace(self):
        a, = self.load("template = 'bodies/test.toml'\nunset = ['atmosphere']\n"
                       "[bodies.atmosphere]\nheight_km = 3\n")
        self.assertEqual(a['atmosphere'], {'height_km': 3})
        self.assertNotIn('unset', a)

    def test_inline_and_empty_table(self):
        self.assertEqual(self.load('radii_km = [1]\n'), [{'name': 'A', 'radii_km': [1]}])
        a, = self.load("template = 'bodies/test.toml'\n[bodies.atmosphere]\n")
        self.assertEqual(a['atmosphere']['gain'], 2)

    def test_bad_references_and_unset(self):
        for value in ("''", '1', "'../test.toml'", "'/test.toml'", "'C:/test.toml'", "'bodies\\test.toml'"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.load(f'template = {value}\n')
        for value in ("'atmosphere'", '[1]', "['orbit']", "['atmosphere', 'atmosphere']"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.load(f"template = 'bodies/test.toml'\nunset = {value}\n")
        with self.assertRaisesRegex(ValueError, 'needs a body template'):
            self.load('unset = []\n')
        with self.assertRaisesRegex(ValueError, "body 'A'.*missing.toml"):
            self.load("template = 'bodies/missing.toml'\n")

    def test_single_level_and_parse_error(self):
        for text in ("template = 'bodies/test.toml'\n", "template = 'bodies/other.toml'\n", 'radii_km = ['):
            self.template.write_text(text, encoding='utf-8')
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, "body 'A'.*bodies/test.toml"):
                self.load("template = 'bodies/test.toml'\n")


if __name__ == '__main__':
    unittest.main()
