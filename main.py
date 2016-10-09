from PyQt4 import QtGui, QtCore
import sys
import design
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""

    def compute_initial_figure(self):
        pass



class ExampleApp(QtGui.QMainWindow, design.Ui_MainWindow):
    def __init__(self, parent=None):
        super(ExampleApp, self).__init__(parent)
        self.setupUi(self)
        self.radioButton_man.toggled.connect(self.disallow_beta_edge)
        self.radioButton_auto.toggled.connect(self.enable_beta_edge)
        self.radioButton_man.setChecked(True)
        self.pushButton.clicked.connect(self.solve)

        self.sol_plot = MyStaticMplCanvas(self.widget_plot, width=5, height=4, dpi=100)

    def enable_beta_edge(self):
        self.lineEdit_beta_to.setEnabled(True)

    def disallow_beta_edge(self):
        self.lineEdit_beta_to.setEnabled(False)

    def error_message(self, param, expected):
        self.label_res.setText("Error: " + param + " must be " + expected + ".")

    def check_params(self):
        for p in self.params:
            if p == "beta_to" and not self.auto_mode:
                continue
            if not self.params[p] or not is_number(self.params[p]):
                self.error_message(p, "a number")
            self.params[p] = float(self.params[p])

    def check_function(self, fun_name, fun):
        if not fun:
            self.error_message(fun_name, "specified")

        t = 0
        try:
            eval(fun)
        except Exception as e:
            self.error_message(fun_name, "valid function of t")

    def check_valid_input(self):
        self.params = {"a": self.lineEdit_a.text(),
                  "b": self.lineEdit_b.text(),
                  "x0": self.lineEdit_x0.text(),
                  "y0": self.lineEdit_y0.text(),
                  "beta_from": self.lineEdit_beta_from.text(),
                  "beta_to": self.lineEdit_beta_to.text()}
        if not self.check_params():
            return False
        self.z = self.lineEdit_z.text()
        self.s = self.lineEdit_s.text()
        if not self.check_function("s", self.s) or not self.check_function("z", self.z):
            return False
        return True

    def solve(self):
        self.label_res.setText("")
        print("Solver called")
        self.auto_mode = self.radioButton_auto.isChecked()
        if not check_valid_input():
            return 
        



    def tabulate(self):
        print("Tabulate called")

    def tabulate_integral(self):
        print("Tabulate integral called")

    def functionFBeta(self):
        print("Threshld tuner called")

    def interpolate(self):
        print("Interpolate called")

    def diff_equation(self):
        print("Diff equation called")

    def lin_sys(self):
        print("Lin sys called")

    def integrate(self):
        print("Integrate called")


def main():
    app = QtGui.QApplication(sys.argv)
    form = ExampleApp()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
