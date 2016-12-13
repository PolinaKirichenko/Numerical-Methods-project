from PyQt4 import QtGui, QtCore
import sys
import new_design
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import numexpr as ne

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

    def solution_plot(self):
        t = np.arange(0.0, 3.0, 0.01)
        s = np.sin(2*np.pi*t)
        self.axes.plot(t, s)
        self.draw()

class ExampleApp(QtGui.QMainWindow, new_design.Ui_MainWindow):
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
                return False
            self.params[p] = float(self.params[p])
            if p == "T" and self.params[p] <= 0:
                self.error_message(p, "a positive number")
                return False
        return True

    def check_function(self, fun_name, fun):
        if not fun:
            self.error_message(fun_name, "specified")
            return False

        t = 0
        try:
            eval(fun)
        except Exception as e:
            self.error_message(fun_name, "valid function of t")
            return False
        return True

    def check_valid_input(self):
        self.params = {"a": self.lineEdit_a.text(),
                  "b": self.lineEdit_b.text(),
                  "x0": self.lineEdit_x0.text(),
                  "y0": self.lineEdit_y0.text(),
                  "T": self.lineEdit_T.text(),
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
        self.grid_size = 10
        if not self.check_valid_input():
            return
        self.target_distr = str(self.params["a"]) + " * t * (" + str(self.params["b"]) + " - t)"
        
        self.tabulate(self.target_distr, "target_distr_tab")
        self.tabulate(self.s, "plan_tab")
        self.tabulate(self.z, "traffic_tab")

        self.tabulate_integral("target_distr_tab")

        self.interpolate("u_tab")
        self.interpolate("plan_tab")
        self.interpolate("traffic_tab")

        #self.diff_equation()

        self.solution_plot("cauchy_sol")

    def solution_plot(self, input_file):
        print("Plotting solution")
        self.sol_plot.solution_plot()

    def tabulate(self, fun, dest_file):
        print("Tabulate called")
        t = np.linspace(0, self.params["T"], self.grid_size)
        np.savetxt(dest_file, (t, ne.evaluate(fun + " + 0 * t")))

    def tabulate_integral(self, input_file):
        print("Tabulate integral called")
        tab_fun = np.genfromtxt(input_file)
        for i in range(tab_fun.shape[0]):
            self.integrate(tab_fun[:, i:])
        np.savetxt("u_tab", np.zeros(self.grid_size))

    def integrate(self, grid):
        print("Integrate")

    def functionFBeta(self):
        print("Threshold tuner called")

    def interpolate(self, input_file):
        print("Interpolate called")
        np.savetxt(input_file[:-4] + "_coef", np.zeros(self.grid_size))

    def diff_equation(self, init_point, integrand):
        ### Euler method ###
        print("Diff equation called")
        np.savetxt("cauchy_sol", np.zeros(self.grid_size))
        return

        grid = np.linspace(0, self.params["T"], self.grid_size)
        sol = np.zeros(self.grid_size)
        for i, t in enumerate(grid):
            if i == 0:
                sol[i] = init_point
            else:
                sol[i] = sol[i - 1] + (grid[i] - grid[i - 1]) * integrand(grid[i - 1], sol[i - 1])
        np.savetxt("cauchy_sol", sol)

    def lin_sys(self):
        print("Lin sys called")


def main():
    app = QtGui.QApplication(sys.argv)
    form = ExampleApp()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
