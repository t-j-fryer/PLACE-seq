"""Supplying the pooling layout, and being told when it is wrong.

The layout is experimental design: which block was picked into which culture
plate, and which plates were pooled into which barcode.  It is known before any
read exists and is never inferred from the data - so the only thing that can go
wrong is transcription, and these tests pin that every transcription error is
named at preflight rather than appearing as unattributable reads three stages on.
"""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from nanopore3.layout import LayoutError, read_layout_csv, write_layout_csv

TABLE = """plate_barcode,culture_plate,library,block
RP05,PlateA,sumo,Block_1
RP05,PlateB,sumo,Block_2
RP05,PlateB,sumo,Block_3
RP08,PlateC,lab,Block_1
"""


class ReadLayoutTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = TemporaryDirectory()
        self.path = Path(self.tmp.name) / "layout.csv"
        self.addCleanup(self.tmp.cleanup)

    def write(self, text: str) -> Path:
        self.path.write_text(text, encoding="utf-8")
        return self.path

    def test_one_row_per_block_and_plate(self) -> None:
        layout = read_layout_csv(self.write(TABLE))
        self.assertEqual(
            layout["pcr_plates"], {"RP05": ("PlateA", "PlateB"), "RP08": ("PlateC",)}
        )
        # Block numbering restarts per library, so the two Block_1s stay apart.
        self.assertEqual(
            layout["blocks"],
            {
                "sumo": {"Block_1": ("PlateA",), "Block_2": ("PlateB",),
                         "Block_3": ("PlateB",)},
                "lab": {"Block_1": ("PlateC",)},
            },
        )

    def test_a_block_split_across_plates_gets_several_rows(self) -> None:
        layout = read_layout_csv(
            self.write(TABLE + "RP05,PlateA,sumo,Block_2\n")
        )
        self.assertEqual(layout["blocks"]["sumo"]["Block_2"], ("PlateB", "PlateA"))

    def test_a_repeated_row_is_not_counted_twice(self) -> None:
        layout = read_layout_csv(self.write(TABLE + "RP05,PlateA,sumo,Block_1\n"))
        self.assertEqual(layout["blocks"]["sumo"]["Block_1"], ("PlateA",))
        self.assertEqual(layout["pcr_plates"]["RP05"], ("PlateA", "PlateB"))

    def test_a_misspelled_column_is_refused_not_read_as_blank(self) -> None:
        with self.assertRaises(LayoutError) as caught:
            read_layout_csv(self.write(TABLE.replace("culture_plate", "culture_plt")))
        self.assertIn("culture_plate", str(caught.exception))

    def test_an_extra_column_is_refused(self) -> None:
        with self.assertRaises(LayoutError) as caught:
            read_layout_csv(self.write(TABLE.replace("block\n", "block,notes\n")))
        self.assertIn("notes", str(caught.exception))

    def test_an_empty_cell_names_its_line(self) -> None:
        with self.assertRaises(LayoutError) as caught:
            read_layout_csv(self.write(TABLE + "RP05,,sumo,Block_9\n"))
        self.assertIn("line 6", str(caught.exception))
        self.assertIn("culture_plate", str(caught.exception))

    def test_blank_lines_and_a_byte_order_mark_are_tolerated(self) -> None:
        # Both are what a spreadsheet export actually produces.
        layout = read_layout_csv(self.write("﻿" + TABLE + "\n\n"))
        self.assertEqual(len(layout["pcr_plates"]), 2)

    def test_an_empty_table_is_an_error_not_an_empty_layout(self) -> None:
        with self.assertRaises(LayoutError):
            read_layout_csv(self.write("plate_barcode,culture_plate,library,block\n"))

    def test_a_missing_file_says_so(self) -> None:
        with self.assertRaises(LayoutError) as caught:
            read_layout_csv(Path(self.tmp.name) / "absent.csv")
        self.assertIn("absent.csv", str(caught.exception))

    def test_clonality_describes_the_barcode_so_it_must_agree(self) -> None:
        text = TABLE.replace(
            "library,block\n", "library,block,clonality\n"
        ).replace("sumo,Block_1\n", "sumo,Block_1,per_block\n", 1)
        text = text.replace("sumo,Block_2\n", "sumo,Block_2,unspecified\n", 1)
        with self.assertRaises(LayoutError) as caught:
            read_layout_csv(self.write(text))
        self.assertIn("RP05", str(caught.exception))

    def test_it_round_trips_through_the_writer(self) -> None:
        layout = read_layout_csv(self.write(TABLE))
        out = Path(self.tmp.name) / "again.csv"
        write_layout_csv(out, layout["pcr_plates"], layout["blocks"])
        self.assertEqual(read_layout_csv(out)["blocks"], layout["blocks"])
        self.assertEqual(read_layout_csv(out)["pcr_plates"], layout["pcr_plates"])


if __name__ == "__main__":
    unittest.main()
