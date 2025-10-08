import logging

def setup_logger(level=logging.INFO):
    logging.basicConfig(level=level,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler("varannototar.log"),
                            logging.StreamHandler()
                        ])