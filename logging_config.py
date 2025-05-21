import logging
import os

def setup_logging():
    log_level = os.getenv("LOG_LEVEL", "INFO").upper()
    log_format = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    log_file = os.getenv("LOG_FILE", "cryoprotectant_analysis.log")

    handlers = [logging.StreamHandler()]
    if os.getenv("LOG_TO_FILE", "1") == "1":
        handlers.append(logging.FileHandler(log_file, mode="a", encoding="utf-8"))

    logging.basicConfig(
        level=getattr(logging, log_level, logging.INFO),
        format=log_format,
        handlers=handlers
    )