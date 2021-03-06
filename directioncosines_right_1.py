# Direction Cosines Right Angular Velocity 1 Program in Python.
# Calculates the Numerical Approximation of the Motion of a 3-D Gyroscopic Pendulum
# Using 1st-Order Runge-Kutta Method (Euler Method).

import math

#
# zyz-Convention Euler Angles Representation of the Proper Orthogonal Linear Transformation
#


def oxx(qyaw, qpitch, qroll):
    f = math.cos(qyaw) * math.cos(qpitch) * math.cos(qroll) - math.sin(qyaw) * math.sin(qroll)
    return f


def oyx(qyaw, qpitch, qroll):
    f = math.sin(qyaw) * math.cos(qpitch) * math.cos(qroll) + math.cos(qyaw) * math.sin(qroll)
    return f


def ozx(qpitch, qroll):
    f = - math.sin(qpitch) * math.cos(qroll)
    return f


def oxy(qyaw, qpitch, qroll):
    f = - math.cos(qyaw) * math.cos(qpitch) * math.sin(qroll) - math.sin(qyaw) * math.cos(qroll)
    return f


def oyy(qyaw, qpitch, qroll):
    f = - math.sin(qyaw) * math.cos(qpitch) * math.sin(qroll) + math.cos(qyaw) * math.cos(qroll)
    return f


def ozy(qpitch, qroll):
    f = math.sin(qpitch) * math.sin(qroll)
    return f


def oxz(qyaw, qpitch):
    f = math.cos(qyaw) * math.sin(qpitch)
    return f


def oyz(qyaw, qpitch):
    f = math.sin(qyaw) * math.sin(qpitch)
    return f


def ozz(qpitch):
    f = math.cos(qpitch)
    return f


#
# zyz-Convention Euler Angles Representation of the Left Angular Velocity
#


def omegax(qpitch, qroll, dqyaw, dqpitch):
    f = - dqyaw * math.sin(qpitch) * math.cos(qroll) + dqpitch * math.sin(qroll)
    return f


def omegay(qpitch, qroll, dqyaw, dqpitch):
    f = dqyaw * math.sin(qpitch) * math.sin(qroll) + dqpitch * math.cos(qroll)
    return f


def omegaz(qpitch, dqyaw, dqroll):
    f = dqyaw * math.cos(qpitch) + dqroll
    return f


def main():
    ns = 700000
    nt = 700
    a = 0.25
    h = 0.1
    l = 1.0
    ipum = [0.0, (a ** 2 + h ** 2 / 3.0) / 4.0 + l ** 2, (a ** 2 + h ** 2 / 3.0) / 4.0 + l ** 2, a ** 2 / 2.0]
    gamma = [0.0, (ipum[2] - ipum[3]) / ipum[1], (ipum[3] - ipum[1]) / ipum[2], (ipum[1] - ipum[2]) / ipum[3]]
    beta = [0.0, l / ipum[1], l / ipum[2], l / ipum[3]]
    bg = [0.0, 0.0, 0.0, -9.8]
    time = [0.0]
    delta = (7.0 - time[0]) / ns
    for i in range(1, ns + 1, 1):
        time.append(time[i-1] + delta)
    yaw = 0.0 * math.pi
    pitch = 11.0 * math.pi / 24.0
    roll = 0.0 * math.pi
    o = [[0.0, 0.0, 0.0, 0.0],
         [0.0, oxx(yaw, pitch, roll), oxy(yaw, pitch, roll), oxz(yaw, pitch)],
         [0.0, oyx(yaw, pitch, roll), oyy(yaw, pitch, roll), oyz(yaw, pitch)],
         [0.0, ozx(pitch, roll), ozy(pitch, roll), ozz(pitch)]]
    bi = [0.0, o[1][1], o[2][1], o[3][1]]
    bj = [0.0, o[1][2], o[2][2], o[3][2]]
    bk = [0.0, o[1][3], o[2][3], o[3][3]]
    bgbi = bg[1] * bi[1] + bg[2] * bi[2] + bg[3] * bi[3]
    bgbj = bg[1] * bj[1] + bg[2] * bj[2] + bg[3] * bj[3]
    bgbk = bg[1] * bk[1] + bg[2] * bk[2] + bg[3] * bk[3]
    dyaw = 0.0 * math.pi
    dpitch = 0.0 * math.pi
    droll = 100.0 * math.pi
    omega = [0.0, omegax(pitch, roll, dyaw, dpitch), omegay(pitch, roll, dyaw, dpitch), omegaz(pitch, dyaw, droll)]
    domega = [0.0,
              gamma[1] * omega[2] * omega[3] - beta[1] * bgbj,
              gamma[2] * omega[3] * omega[1] + beta[2] * bgbi,
              gamma[3] * omega[1] * omega[2]]
    oj = o.copy()
    omegaj = omega.copy()
    ot = [[[ 0.0 for j in range(4)] for i in range(4)] for k in range(nt + 1)]
    ot[0][1][1] = o[1][1]
    ot[0][2][1] = o[2][1]
    ot[0][3][1] = o[3][1]
    ot[0][1][2] = o[1][2]
    ot[0][2][2] = o[2][2]
    ot[0][3][2] = o[3][2]
    ot[0][1][3] = o[1][3]
    ot[0][2][3] = o[2][3]
    ot[0][3][3] = o[3][3]
    omegat = [[ 0.0 for i in range(4)] for k in range(nt + 1)]
    omegat[0][1] = omega[1]
    omegat[0][2] = omega[2]
    omegat[0][3] = omega[3]
    bit = [[ 0.0 for i in range(4)] for k in range(nt + 1)]
    bit[0][1] = bi[1]
    bit[0][2] = bi[2]
    bit[0][3] = bi[3]
    bjt = [[ 0.0 for i in range(4)] for k in range(nt + 1)]
    bjt[0][1] = bj[1]
    bjt[0][2] = bj[2]
    bjt[0][3] = bj[3]
    bkt = [[ 0.0 for i in range(4)] for k in range(nt + 1)]
    bkt[0][1] = bk[1]
    bkt[0][2] = bk[2]
    bkt[0][3] = bk[3]
    for k in range(1, nt + 1, 1):
        for j in range(1, ns // nt + 1, 1):
            nj = ns // nt * (k - 1) + j - 1
            deltat = time[nj + 1] - time[nj]
            fo1 = [[0.0, 0.0, 0.0, 0.0],
                   [0.0, o[1][2] * omega[3] - o[1][3] * omega[2],
                    o[1][3] * omega[1] - o[1][1] * omega[3],
                    o[1][1] * omega[2] - o[1][2] * omega[1]],
                   [0.0, o[2][2] * omega[3] - o[2][3] * omega[2],
                    o[2][3] * omega[1] - o[2][1] * omega[3],
                    o[2][1] * omega[2] - o[2][2] * omega[1]],
                   [0.0, o[3][2] * omega[3] - o[3][3] * omega[2],
                    o[3][3] * omega[1] - o[3][1] * omega[3],
                    o[3][1] * omega[2] - o[3][2] * omega[1]]]
            fomega1 = domega.copy()
            o[1][1] = oj[1][1] + deltat * fo1[1][1]
            o[2][1] = oj[2][1] + deltat * fo1[2][1]
            o[3][1] = oj[3][1] + deltat * fo1[3][1]
            o[1][2] = oj[1][2] + deltat * fo1[1][2]
            o[2][2] = oj[2][2] + deltat * fo1[2][2]
            o[3][2] = oj[3][2] + deltat * fo1[3][2]
            o[1][3] = oj[1][3] + deltat * fo1[1][3]
            o[2][3] = oj[2][3] + deltat * fo1[2][3]
            o[3][3] = oj[3][3] + deltat * fo1[3][3]
            omega[1] = omegaj[1] + deltat * fomega1[1]
            omega[2] = omegaj[2] + deltat * fomega1[2]
            omega[3] = omegaj[3] + deltat * fomega1[3]
            bi = [0.0, o[1][1], o[2][1], o[3][1]]
            bj = [0.0, o[1][2], o[2][2], o[3][2]]
            bk = [0.0, o[1][3], o[2][3], o[3][3]]
            bgbi = bg[1] * bi[1] + bg[2] * bi[2] + bg[3] * bi[3]
            bgbj = bg[1] * bj[1] + bg[2] * bj[2] + bg[3] * bj[3]
            bgbk = bg[1] * bk[1] + bg[2] * bk[2] + bg[3] * bk[3]
            domega = [0.0,
                      gamma[1] * omega[2] * omega[3] - beta[1] * bgbj,
                      gamma[2] * omega[3] * omega[1] + beta[2] * bgbi,
                      gamma[3] * omega[1] * omega[2]]
            oj = o.copy()
            omegaj = omega.copy()
        ot[k][1][1] = o[1][1]
        ot[k][2][1] = o[2][1]
        ot[k][3][1] = o[3][1]
        ot[k][1][2] = o[1][2]
        ot[k][2][2] = o[2][2]
        ot[k][3][2] = o[3][2]
        ot[k][1][3] = o[1][3]
        ot[k][2][3] = o[2][3]
        ot[k][3][3] = o[3][3]
        omegat[k][1] = omega[1]
        omegat[k][2] = omega[2]
        omegat[k][3] = omega[3]
        bit[k][1] = bi[1]
        bit[k][2] = bi[2]
        bit[k][3] = bi[3]
        bjt[k][1] = bj[1]
        bjt[k][2] = bj[2]
        bjt[k][3] = bj[3]
        bkt[k][1] = bk[1]
        bkt[k][2] = bk[2]
        bkt[k][3] = bk[3]
        print("%7u%3c%- 22.15e" % (nj + 1, ' ', time[nj + 1]))
    unit1 = open("dc_r_1_primary.out", "w")
    unit2 = open("dc_r_1_secondary.out", "w")
    unit1.write("%9c%4s%20c%6s%19c%6s%19c%6s%21c%3s%22c%3s%22c%3s%22c%3s%22c%3s%22c%3s%22c%3s%22c%3s%22c%3s\n" %
                (' ', "time", ' ', "omegax", ' ', "omegay", ' ', "omegaz", ' ', "bix", ' ', "biy", ' ', "biz",
                 ' ', "bjx", ' ', "bjy", ' ', "bjz", ' ', "bkx", ' ', "bky", ' ', "bkz"))
    unit2.write("%9c%4s%21c%5s%20c%5s%20c%5s%20c%4s%21c%4s%21c%4s%21c%4s\n" %
                (' ', "time", ' ', "bimag", ' ', "bjmag", ' ', "bkmag", ' ', "bibj", ' ', "bkbi", ' ', "bjbk", ' ',
                 "epum"))
    for k in range(0, nt + 1, 1):
        bimag = math.sqrt(bit[k][1] ** 2 + bit[k][2] ** 2 + bit[k][3] ** 2)
        bjmag = math.sqrt(bjt[k][1] ** 2 + bjt[k][2] ** 2 + bjt[k][3] ** 2)
        bkmag = math.sqrt(bkt[k][1] ** 2 + bkt[k][2] ** 2 + bkt[k][3] ** 2)
        bibj = bit[k][1] * bjt[k][1] + bit[k][2] * bjt[k][2] + bit[k][3] * bjt[k][3]
        bkbi = bkt[k][1] * bit[k][1] + bkt[k][2] * bit[k][2] + bkt[k][3] * bit[k][3]
        bjbk = bjt[k][1] * bkt[k][1] + bjt[k][2] * bkt[k][2] + bjt[k][3] * bkt[k][3]
        bgbk = bg[1] * bkt[k][1] + bg[2] * bkt[k][2] + bg[3] * bkt[k][3]
        epum = (ipum[1] * omegat[k][1] ** 2 + ipum[2] * omegat[k][2] ** 2 + ipum[3] * omegat[k][3] ** 2) / 2.0\
               - l * bgbk
        unit1.write(
         "%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c"
         "%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e\n"
         % (time[ns // nt * k], ' ',
            omegat[k][1], ' ', omegat[k][2], ' ', omegat[k][3], ' ',
            bit[k][1], ' ', bit[k][2], ' ', bit[k][3], ' ',
            bjt[k][1], ' ', bjt[k][2], ' ', bjt[k][3], ' ',
            bkt[k][1], ' ', bkt[k][2], ' ', bkt[k][3]))
        unit2.write(
         "%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e%3c%- 22.15e\n"
         % (time[ns // nt * k], ' ', bimag, ' ', bjmag, ' ', bkmag, ' ', bibj, ' ', bkbi, ' ', bjbk, ' ', epum))
    unit1.close()
    unit2.close()


if __name__ == '__main__':
    main()
    exit(0)
