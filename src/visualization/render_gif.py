#!/usr/bin/env python3
"""Render a time-resolved VTK electron scene and assemble its frames as GIF."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gif", type=Path, required=True)
    parser.add_argument("--frames", type=int, default=80)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--width", type=int, default=960)
    parser.add_argument("--height", type=int, default=600)
    parser.add_argument("--visualizer", type=Path, default=
                        Path(__file__).resolve().parents[2] / "build/aurora_vtk_visualizer")
    args, forwarded = parser.parse_known_args()
    if args.gif.suffix.lower() != ".gif" or args.frames < 2 or args.frames > 600:
        parser.error("--gif must end in .gif and --frames must be in [2,600]")
    if args.fps < 1 or args.fps > 50:
        parser.error("--fps must be in [1,50]")
    if any(option in forwarded for option in ("--output", "--frames-dir",
                                                "--frames", "--width", "--height")):
        parser.error("Pass frame size/count to this script, not to the visualizer")
    try:
        from PIL import Image, ImageDraw
    except ImportError as error:
        raise SystemExit("Pillow is required: python3 -m pip install Pillow") from error

    destination = args.gif.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    preview = destination.with_suffix(".png")
    paired = "--compare-anomaly" in forwarded
    with tempfile.TemporaryDirectory(prefix="aurora-gif-frames-") as temporary:
        command = [str(args.visualizer), "--output", str(preview),
                   "--frames-dir", temporary, "--frames", str(args.frames),
                   "--width", str(args.width), "--height", str(args.height),
                   *forwarded]
        subprocess.run(command, check=True)
        if paired:
            base_frames = sorted((Path(temporary) / "base").glob("frame_*.png"))
            anomaly_frames = sorted((Path(temporary) / "anomaly").glob("frame_*.png"))
            if len(base_frames) != args.frames or len(anomaly_frames) != args.frames:
                raise RuntimeError("Paired renderer did not produce all frames")
            frame_pairs = zip(base_frames, anomaly_frames)
        else:
            frame_paths = sorted(Path(temporary).glob("frame_*.png"))
            if len(frame_paths) != args.frames:
                raise RuntimeError(f"Expected {args.frames} frames, got {len(frame_paths)}")
            frame_pairs = ((path, None) for path in frame_paths)
        images = []
        for base_path, anomaly_path in frame_pairs:
            with Image.open(base_path) as source:
                frame = source.convert("RGB")
            if anomaly_path is not None:
                with Image.open(anomaly_path) as source:
                    anomaly = source.convert("RGB")
                combined = Image.new("RGB", (frame.width + anomaly.width,
                                             max(frame.height, anomaly.height) + 36),
                                     (8, 16, 33))
                combined.paste(frame, (0, 36))
                combined.paste(anomaly, (frame.width, 36))
                labels = ImageDraw.Draw(combined)
                labels.text((16, 10), "BASE DIPOLE", fill=(230, 240, 255))
                labels.text((frame.width + 16, 10), "WITH ANOMALY",
                            fill=(230, 240, 255))
                frame = combined
            images.append(frame.quantize(colors=128))
        if paired:
            images[-1].convert("RGB").save(preview)
        temporary_gif = destination.with_name(destination.name + ".partial")
        try:
            images[0].save(temporary_gif, format="GIF", save_all=True,
                           append_images=images[1:], duration=round(1000 / args.fps),
                           loop=0, disposal=2, optimize=False)
            os.replace(temporary_gif, destination)
        finally:
            temporary_gif.unlink(missing_ok=True)
    print(f"Saved {destination} and final-frame preview {preview}")


if __name__ == "__main__":
    main()
