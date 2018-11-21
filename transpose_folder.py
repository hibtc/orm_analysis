import os
import sys

src_dir = sys.argv[1]
dst_dir = sys.argv[2]
select = sys.argv[3]

files = [
    (d, f)
    for d in os.listdir(src_dir)
    if os.path.isdir(os.path.join(src_dir, d))
    for f in os.listdir(os.path.join(src_dir, d))
    if f.endswith(select)
]

for d in {os.path.join(dst_dir, f) for _, f in files}:
    os.makedirs(d, exist_ok=True)

for d, f in files:
    os.link(
        os.path.join(src_dir, d, f),
        os.path.join(dst_dir, f, d) + os.path.splitext(f)[1])
