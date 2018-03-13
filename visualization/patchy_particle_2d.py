import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def plot_Janus(ax,
               xc,
               yc,
               radius,
               theta,
               delta=None,
               color=["r", "b"],
               resolution=100):
    if delta is None:
        delta = np.pi / 2
        half_circle = True
    else:
        half_circle = False
    n1 = int(resolution * delta / np.pi)
    theta1 = np.linspace(theta - delta, theta + delta, n1)
    xy1 = np.zeros((n1, 2))
    xy1[:, 0] = radius * np.cos(theta1) + xc
    xy1[:, 1] = radius * np.sin(theta1) + yc
    ax.add_patch(patches.Polygon(xy1, fc=color[0]))

    if half_circle:
        xy2 = np.zeros_like(xy1)
        xy2[:, 0] = 2 * xc - xy1[:, 0]
        xy2[:, 1] = 2 * yc - xy1[:, 1]
    else:
        n2 = resolution - n1
        theta2 = np.linspace(theta + delta, theta - delta + np.pi * 2, n2)
        print("n2 =", n2, theta2.size)
        xy2 = np.zeros((n2, 2))
        xy2[:, 0] = radius * np.cos(theta2) + xc
        xy2[:, 1] = radius * np.sin(theta2) + yc
    ax.add_patch(patches.Polygon(xy2, fc=color[1]))


def plot_triblock(ax,
                  xc,
                  yc,
                  radius,
                  theta,
                  delta,
                  color=["b", "y", "g"],
                  resolution=100):
    if isinstance(delta, list):
        delta_head, delta_tail = delta
        n1 = int(resolution * delta_head / np.pi)
        theta1 = np.linspace(theta - delta_head, theta + delta_head, n1)
        xy1 = np.zeros((n1, 2))
        xy1[:, 0] = radius * np.cos(theta1) + xc
        xy1[:, 1] = radius * np.sin(theta1) + yc
        n3 = int(resolution * delta_tail / np.pi)
        theta3 = np.linspace(theta - delta_tail, theta + delta_tail,
                             n3) + np.pi
        xy3 = np.zeros((n3, 2))
        xy3[:, 0] = np.cos(theta3) + xc
        xy3[:, 1] = np.sin(theta3) + yc
        n2 = resolution - n1 - n3
        n2_half = int(n2 / 2)
        theta2 = np.zeros(n2)
        theta2[:n2_half] = np.linspace(theta - np.pi + delta_tail,
                                       theta - delta_head, n2_half)
        theta2[n2_half:] = np.linspace(
            theta + delta_head, theta + np.pi - delta_tail, n2 - n2_half)
        xy2 = np.zeros((n2, 2))
        xy2[:, 0] = radius * np.cos(theta2) + xc
        xy2[:, 1] = radius * np.sin(theta2) + yc
    else:
        n1 = int(resolution * delta / np.pi)
        theta1 = np.linspace(theta - delta, theta + delta, n1)
        xy1 = np.zeros((n1, 2))
        xy1[:, 0] = radius * np.cos(theta1) + xc
        xy1[:, 1] = radius * np.sin(theta1) + yc
        xy3 = np.zeros_like(xy1)
        xy3[:, 0] = 2 * xc - xy1[:, 0]
        xy3[:, 1] = 2 * yc - xy1[:, 1]
        half_n2 = int((resolution - 2 * n1) / 2)
        xy2 = np.zeros((half_n2 * 2, 2))
        theta2 = np.linspace(theta - np.pi + delta, theta - delta, half_n2)
        xy2[:half_n2, 0] = radius * np.cos(theta2) + xc
        xy2[:half_n2, 1] = radius * np.sin(theta2) + yc
        xy2[half_n2:, 0] = 2 * xc - xy2[:half_n2, 0]
        xy2[half_n2:, 1] = 2 * yc - xy2[:half_n2, 1]
    ax.add_patch(patches.Polygon(xy1, fc=color[0]))
    ax.add_patch(patches.Polygon(xy2, fc=color[1]))
    if isinstance(delta, list):
        ax.add_patch(patches.Polygon(xy3, fc=color[2]))
    else:
        ax.add_patch(patches.Polygon(xy3, fc=color[0]))


class MgCu2_110:
    def __init__(self, r_l=np.sqrt(3 / 2)):
        self.r_s = 1
        self.r_l = r_l
        self.get_unit_cell()

    def get_unit_cell(self):
        self.vec_a = np.array([self.r_s * 4, 0])
        self.p_s = np.array([[0, 0], [2 * self.r_s, 0]])

        xl_0 = self.r_s
        yl_0 = np.sqrt((self.r_s + self.r_l)**2 - self.r_s**2)
        xl_1 = 3 * self.r_s
        yl_1 = np.sqrt(self.r_l**2 * 4 - self.r_s**2 * 4) + yl_0
        self.vec_b = np.array([self.r_s * 4, yl_0 + yl_1])
        self.p_l = np.array([[xl_0, yl_0], [xl_1, yl_1]])
        self.theta_s = np.arccos(self.r_s / (self.r_s + self.r_l))
        print(self.theta_s / np.pi * 180)
        self.theta_l = np.pi / 2

    def show(self, na=1, nb=1, delta_s=None, delta_l=None):
        ax = plt.subplot(111)
        for ia in range(na + 1):
            for ib in range(nb + 1):
                vec = self.vec_a * ia + self.vec_b * ib
                if delta_s is None:
                    ax.add_patch(patches.Circle(self.p_s[0] + vec, self.r_s))
                    ax.add_patch(patches.Circle(self.p_s[1] + vec, self.r_s))
                else:
                    plot_triblock(ax, self.p_s[0][0] + vec[0],
                                  self.p_s[0][1] + vec[1], self.r_s,
                                  self.theta_s, delta_s)
                    plot_triblock(ax, self.p_s[1][0] + vec[0],
                                  self.p_s[1][1] + vec[1], self.r_s,
                                  np.pi - self.theta_s, delta_s)
                if delta_l is None:
                    ax.add_patch(
                        patches.Circle(self.p_l[0] + vec, self.r_l, fc="r"))
                    ax.add_patch(
                        patches.Circle(self.p_l[1] + vec, self.r_l, fc="r"))
                else:
                    plot_Janus(ax, self.p_l[0][0] + vec[0],
                               self.p_l[0][1] + vec[1], self.r_l, self.theta_l,
                               delta_l)
                    plot_Janus(ax, self.p_l[1][0] + vec[0],
                               self.p_l[1][1] + vec[1], self.r_l,
                               -self.theta_l, delta_l)
        ax.axis("equal")
        plt.show()
        plt.close()


if __name__ == "__main__":
    # ax = plt.subplot(111)
    # plot_Janus(ax, 0.5, 0.5, 1, np.pi / 4)
    # plot_triblock(ax, 3, 3, np.sqrt(2 / 3), -np.pi / 4, np.pi / 3)
    # ax.axis("equal")
    # plt.show()
    # plt.close()
    plane = MgCu2_110()
    plane.show(3, 2, np.pi / 4, np.pi/2)
    # plane.show(3, 2)
