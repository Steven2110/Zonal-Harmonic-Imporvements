import logging

class Logger:
    def __init__(self, log_file):
        """
        Initialize a logger instance with a specified log file.
        
        :param log_file: Path to the log file.
        """
        self.logger = logging.getLogger(log_file)
        self.logger.setLevel(logging.INFO)

        # Prevent duplicate log entries
        if not self.logger.handlers:
            # Create a file handler
            file_handler = logging.FileHandler(log_file, mode='w')
            file_handler.setLevel(logging.INFO)

            # Create a logging format
            formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
            file_handler.setFormatter(formatter)

            # Add the file handler to the logger
            self.logger.addHandler(file_handler)

    def log_info(self, message):
        """Logs an INFO level message."""
        self.logger.info(message)

    def log_error(self, message):
        """Logs an ERROR level message."""
        self.logger.error(message)

    def log_warning(self, message):
        """Logs a WARNING level message."""
        self.logger.warning(message)

    def log_debug(self, message):
        """Logs a DEBUG level message."""
        self.logger.debug(message)
