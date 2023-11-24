import os
import sys

os.environ['R_HOME'] = 'C:/Program Files/R/R-4.3.2'

from PyQt5.QtCore import QAbstractListModel, Qt
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout,
                             QTabWidget)
from src.window.preprocessing_window import PreProcessingWindow
from src.window.repurposing_window import RepurposingWindow


class MyListModel(QAbstractListModel):
    def __init__(self, data, parent=None):
        super(MyListModel, self).__init__(parent)
        self._data = data

    def rowCount(self, parent=None, *args, **kwargs):
        return len(self._data)

    def data(self, index, role=None):
        if role == Qt.DisplayRole:
            return self._data[index.row()]


class MyMainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.init_ui()

    def init_ui(self):
        # Create a menu bar
        # menu_bar = self.create_menu_bar()

        # Create a QTabWidget to switch between tabs
        tab_widget = QTabWidget()

        # Create Tab 1 content
        tab1_content = RepurposingWindow()
        tab_widget.addTab(tab1_content, 'Recupero farmaci')  # Placeholder for the first tab

        # Create Tab 2 content by replicating Tab 1 content
        tab2_content = PreProcessingWindow()
        tab_widget.addTab(tab2_content, 'Preprocessamento nuove linee')  # Placeholder for the second tab

        # Connect the button click event to a function

        # Connect the button click event in Tab 2 to a function

        # Create a vertical layout for the main window
        main_window_layout = QVBoxLayout()
        # main_window_layout.addWidget(menu_bar)
        main_window_layout.addWidget(tab_widget)

        # Set the main layout for the main window
        self.setLayout(main_window_layout)

        # Set window properties
        self.setWindowTitle('TSR for drug repurposing')
        self.setGeometry(100, 100, 800, 400)  # (x, y, width, height)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MyMainWindow()
    window.show()
    sys.exit(app.exec_())
