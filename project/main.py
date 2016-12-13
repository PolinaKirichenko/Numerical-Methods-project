from PyQt4 import QtGui, QtCore
import sys
import fin_design
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import numexpr as ne
import num_integral
import spline
import runge
from math import *

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

    def solution_plot(self, x, y):
        self.axes.plot(x, y)
        self.draw()

    def mult_solution_plot(self, xs, ys):
        self.axes.plot(xs[0], ys[0], color='blue')
        self.axes.hold(True)
        for i in range(1, len(xs)):
            self.axes.plot(xs[i], ys[i], color='blue')
        self.axes.hold(False)
        self.draw()

class ExampleApp(QtGui.QMainWindow, fin_design.Ui_MainWindow):
    def __init__(self, parent=None):
        super(ExampleApp, self).__init__(parent)
        self.setupUi(self)
        self.radioButton_man.toggled.connect(self.disallow_beta_edge)
        self.radioButton_auto.toggled.connect(self.enable_beta_edge)
        self.radioButton_man.setChecked(True)
        self.pushButton.clicked.connect(self.solve)
        self.pushButton_plot.setEnabled(False)
        self.pushButton_plot.clicked.connect(self.solution_plot)

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


    def C1(self, x, y):
        t_grid = np.linspace(0, self.params["T"], self.grid_size)
        y_grid = np.vectorize(y.at)(t_grid)
        dx_grid = np.vectorize(x.grad)(t_grid)
        int_grid = []
        f = lambda w: w * np.vectorize(self.spl_ro.at)(w)
        for _y in y_grid:
            w_grid = np.linspace(_y, 1, self.grid_size)
            int_grid.append(
                num_integral.integrate_trapez(np.stack((w_grid, f(w_grid))))
            )

        agg = num_integral.integrate_trapez(np.stack((t_grid, dx_grid * np.array(int_grid))))
        return 1 - agg / (x.at(self.params["T"]) - self.params["x0"])


    def C2(self, x, S):
        T = self.params["T"]
        return abs(x.at(T) - S.at(T)) / S.at(T)


    def getXY(self, beta, x0, y0):
        def f(t, x, y):
            tmp1 = self.spl_Z.grad(t)
            tmp2 = self.spl_I.at(y)
            if tmp1 == None or tmp2 == None:
                self.label_res.setText("Bad initial conditions")
                raise
            return tmp1 * tmp2

        def g(t, x, y):
            return self.functionFBeta(beta, t, x)

        tempx, tempy = runge.runge_kutta(self.params["T"], self.grid_size, x0, y0, f, g)

        t_grid = np.linspace(0, self.params["T"], num=self.grid_size)

        sol_x = spline.CubicSpline(t_grid, tempx)
        sol_y = spline.CubicSpline(t_grid, tempy)

        return sol_x, sol_y        


    def solve(self):
        self.grid_size = 100
        self.betaN = 5
        self.init_eps = 0.1
        self.initN = 3
        A = 1
        B = 10

        self.label_res.setText("")
        print("Solver called")
        self.auto_mode = self.radioButton_auto.isChecked()
        if not self.check_valid_input():
            return
        self.target_distr = str(self.params["a"]) + " * t * (" + str(self.params["b"]) + " - t)"
        
        self.tabulate(self.target_distr, "target_distr_tab")
        self.tabulate(self.s, "plan_tab")
        self.tabulate(self.z, "traffic_tab")
        self.tabulate_integral("target_distr_tab", "integral_tab")  

        self.spl_ro = self.interpolate("target_distr_tab")
        self.spl_S = self.interpolate("plan_tab")
        self.spl_Z = self.interpolate("traffic_tab")
        self.spl_I = self.interpolate("integral_tab")

        if self.auto_mode:
            beta_search = np.linspace(self.params["beta_from"], self.params["beta_to"], self.betaN)
        else:
            beta_search = np.array([self.params["beta_from"]])

        fun_val = []
        for beta in beta_search:
            sol_x, sol_y = self.getXY(beta, self.params["x0"], self.params["y0"])
            fun_val.append(A * self.C1(sol_x, sol_y) + B * self.C2(sol_x, self.spl_S))

        self.beta_opt = beta_search[np.argmin(np.array(fun_val))]
        self.sol_x, self.sol_y = self.getXY(self.beta_opt, self.params["x0"], self.params["y0"])

        self.label_res.setText("Beta:" + str(self.beta_opt) +\
                             "\nC1:" + str(self.C1(self.sol_x, self.sol_y)) +\
                             "\nC2:" + str(self.C2(self.sol_x, self.spl_S)))
        
        self.pushButton_plot.setEnabled(True)


    def solution_plot(self):
        cur = self.comboPlot.currentText()
        if cur == "ρ(w)":
            wgrid = np.linspace(0, 1, self.grid_size)
            self.sol_plot.solution_plot(wgrid, np.vectorize(self.spl_ro.at)(wgrid))
            return

        tgrid = np.linspace(0, self.params["T"], self.grid_size)
        d0 = {"z(t)": np.vectorize(self.spl_Z.at),
              "S(t)": np.vectorize(self.spl_S.at),
        }
        if cur in ["z(t)", "S(t)"]:
            self.sol_plot.solution_plot(tgrid, d0[cur](tgrid))
            return

        if not self.auto_mode:
            def fun_diff(t):
                return self.sol_x.at(t) - self.spl_S.at(t)

            if cur == "S(x)":
                self.sol_plot.solution_plot(
                    np.vectorize(self.sol_x.at)(tgrid), np.vectorize(self.spl_S.at)(tgrid))
                return
            d1 = {"y(t)": np.vectorize(self.sol_y.at),
                  "x(t)": np.vectorize(self.sol_x.at),
                  "x(t) - S(t)": np.vectorize(fun_diff),
            }
            self.sol_plot.solution_plot(tgrid, d1[cur](tgrid))
            return 

        else:
            plot_args = []
            plots_vals = []

            x0 = self.params["x0"]; y0 = self.params["y0"]
            Xs = np.linspace(x0, x0 + self.init_eps, self.initN)
            Ys = np.linspace(y0, y0 + self.init_eps, self.initN)
            X, Y = np.meshgrid(Xs, Ys)
            
            for (xi, yi) in np.vstack((X.flatten(), Y.flatten())).T:
                sol_x, sol_y = self.getXY(self.beta_opt, xi, yi)
                if cur == "S(x)":
                    plots_vals.append(np.vectorize(self.spl_S.at)(tgrid))
                    plot_args.append(np.vectorize(sol_x.at)(tgrid))
                    continue
                def fun_diff(t):
                    return sol_x.at(t) - self.spl_S.at(t)
                d1 = {"y(t)": np.vectorize(sol_y.at),
                      "x(t)": np.vectorize(sol_x.at),
                      "x(t) - S(t)": np.vectorize(fun_diff),
                }
                plots_vals.append(d1[cur](tgrid))
                plot_args.append(tgrid)
            self.sol_plot.mult_solution_plot(plot_args, plots_vals)



    def tabulate(self, fun, dest_file):
        #print("Tabulate called")
        t = np.linspace(0, self.params["T"], num=self.grid_size)
        np.savetxt(dest_file, (t, ne.evaluate(fun + " + 0 * t")))

    def tabulate_integral(self, input_file, output_file):
        #print("Tabulate integral called")
        tab_fun = np.genfromtxt(input_file)
        args = tab_fun[0]
        vals = np.zeros(args.size)
        for i in range(tab_fun.shape[0]):
            vals[i] = self.integrate(tab_fun[:, i:])
        np.savetxt(output_file, (args, vals))

    def integrate(self, tab_fun):
        #print("Integrate")
        return num_integral.integrate_trapez(tab_fun)

    def functionFBeta(self, beta, t, x):
        #print("Threshold tuner called")
        tmp = self.spl_S.at(t)
        if tmp == None:
            self.label_res.setText("Bad initial conditions")
            raise
        return beta * (tmp - x)

    def interpolate(self, input_file):
        #print("Interpolate called")
        tab_fun = np.genfromtxt(input_file)
        return spline.CubicSpline(tab_fun[0], tab_fun[1])

    def diff_equation(self, init_point, integrand):
        ### Euler method ###
        #print("Diff equation called")
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


def main():
    app = QtGui.QApplication(sys.argv)
    form = ExampleApp()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
