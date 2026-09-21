import argparse
import os
import time
from datetime import datetime

import cv2
from picamera2 import Picamera2


CAMERA_SIZE = (820, 616)
CAMERA_FPS = 50
EXPOSURE_US = 10000
GAIN_CAM = 12.0
WARMUP_SECONDS = 1.0
VIDEO_SECONDS = 10.0

def main():
    parser = argparse.ArgumentParser(
        description="Prend une photo puis une vidéo avec la caméra."
    )
    parser.add_argument(
        "--outdir", default="data/camera_capture", help="Dossier de sortie."
    )
    parser.add_argument(
        "--duration", type=float, default=VIDEO_SECONDS,
        help="Durée de la vidéo en secondes."
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    photo_path = os.path.join(args.outdir, f"photo_{stamp}.jpg")
    video_path = os.path.join(args.outdir, f"video_{stamp}.mp4")

    camera = Picamera2()
    configuration = camera.create_video_configuration(
        main={"size": CAMERA_SIZE, "format": "RGB888"},
        controls={
            "FrameDurationLimits": (
                int(1_000_000 / CAMERA_FPS), int(1_000_000 / CAMERA_FPS)
            ),
            "ExposureTime": EXPOSURE_US,
            "AnalogueGain": GAIN_CAM,
        },
    )
    camera.configure(configuration)

    writer = None
    try:
        print("[cam] démarrage...")
        camera.start()
        time.sleep(WARMUP_SECONDS)

        frame = camera.capture_array()
        if not cv2.imwrite(photo_path, frame):
            raise RuntimeError(f"Impossible d'enregistrer {photo_path}")
        print(f"[photo] {photo_path}")

        height, width = frame.shape[:2]
        writer = cv2.VideoWriter(
            video_path,
            cv2.VideoWriter_fourcc(*"mp4v"),
            CAMERA_FPS,
            (width, height),
        )
        if not writer.isOpened():
            raise RuntimeError("Impossible d'ouvrir le fichier vidéo.")

        print(f"[vidéo] enregistrement pendant {args.duration:.1f} s...")
        start = time.perf_counter()
        frames = 0
        while time.perf_counter() - start < args.duration:
            writer.write(camera.capture_array())
            frames += 1

        print(f"[vidéo] {video_path} ({frames} images)")
    finally:
        if writer is not None:
            writer.release()
        camera.stop()
        camera.close()
        print("[cam] arrêtée")


if __name__ == "__main__":
    main()
