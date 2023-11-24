from PyQt5.QtCore import QAbstractTableModel, Qt
from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QTableView, QMessageBox, QPushButton, \
    QFileDialog, QHeaderView


class DataFrameViewer(QMainWindow):
    def __init__(self, dataframe, data):
        super().__init__()
        self.data = data
        self.init_ui(dataframe)

    def init_ui(self, dataframe):
        self.setWindowTitle("Tessuto ::" + self.data['selected_tissue'] +
                            " || CCLE::" + self.data["selected_element"] +
                            " || Method::" + self.data["selected_method"] +
                            " || Timeline::" + self.data["selected_timeline"])
        self.setGeometry(100, 100, 1224,800)

        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout()

        # Create a QTableView to display the DataFrame
        table_view = QTableView(self)
        model = PandasModel(dataframe)
        table_view.setModel(model)

        # Set the resize mode for the horizontal header
        table_view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        # Enable sorting for the table view
        table_view.sortByColumn(1,Qt.AscendingOrder)
        table_view.setSortingEnabled(True)

        layout.addWidget(table_view)
        central_widget.setLayout(layout)

        # Create a "Save as CSV" button
        save_button = QPushButton('Save as CSV', self)
        save_button.clicked.connect(self.save_as_csv)
        layout.addWidget(save_button)

    def save_as_csv(self):
        # Get the DataFrame from the model
        model = self.centralWidget().layout().itemAt(0).widget().model()
        original_df = model._data
        self.save_df_as_csv(original_df)

    def save_df_as_csv(self, dataframe):
        # Show a QMessageBox with the DataFrame content
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle('DataFrame Content')
        msg_box.setText("Salvataggio Effettuato!")

        save_file_dialog = QFileDialog(self)
        save_file_dialog.setDefaultSuffix("csv")

        save_file_dialog.setNameFilter("CSV Files (*.csv)")

        if save_file_dialog.exec_():
            save_file_path = save_file_dialog.selectedFiles()[0]
            dataframe.to_csv(save_file_path, index=False)
            msg_box.exec_()
            self.close()


class PandasModel(QAbstractTableModel):
    def __init__(self, data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return len(self._data.index)

    def columnCount(self, parent=None):
        return self._data.columns.size

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return str(self._data.columns[section])
        return None

    def sort(self, column, order):
        self.layoutAboutToBeChanged.emit()
        if order == Qt.AscendingOrder:
            self._data.sort_values(by=self._data.columns[column], inplace=True, ascending=True)
        elif order == Qt.DescendingOrder:
            self._data.sort_values(by=self._data.columns[column], inplace=True, ascending=False)
        self.layoutChanged.emit()