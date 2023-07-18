from PyQt5 import uic, QtWidgets
import sys

path = "/Users/akhiltayal/Library/CloudStorage/GoogleDrive-akhil.tayal.bnl@gmail.com/My Drive/Demo/"

class Shell(*uic.loadUiType(path + 'shell.ui')):
    def __init__(self):
        super().__init__()
        self.setupUi(self)


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    window = Shell()
    window.show()
    sys.exit(app.exec_())