#!/usr/bin/env python3
"""Preserve Figure1 physical panel sizes when included in the paper."""
import sys
from pathlib import Path
from pypdf import PdfReader,PdfWriter
from pypdf.generic import RectangleObject
path=Path(sys.argv[1])
reader=PdfReader(path)
if len(reader.pages)!=1:
    raise ValueError('Expected one Figure1 page')
writer=PdfWriter();writer.clone_document_from_reader(reader)
page=writer.pages[0]
page.mediabox=RectangleObject([0,0,532*2/3,216*2/(2+.17)])
page.cropbox=page.mediabox
with path.open('wb') as stream:
    writer.write(stream)
