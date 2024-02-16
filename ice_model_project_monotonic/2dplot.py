import numpy as np
import matplotlib.pyplot as pl

a = open('t_N.dat', 'r')
b = open('ub_lb_t.dat', 'r')
c = open('T_profile.dat', 'r')
d = open('Tu_Tl.dat', 'r')
f = open('t_xf3.dat', 'r')
g = open('t_xf2.dat', 'r')
h = open('t1.dat', 'r')
l = open('t_dry.dat', 'r')

a1 = a.readlines()
b1 = b.readlines()
c1 = c.readlines()
d1 = d.readlines()
f1 = f.readlines()
g1 = g.readlines()
h1 = h.readlines()
l1 = l.readlines()

del a
del b
del c
del d
del f
del g
del h
del l

time3 = []
time2 = []
time1 = []
N_lens = []
Tl = []
T_ub = []
x_ub = []
x_lb = []
T = []
xf2 = []
xf3 = []
x_dry = []

for k in range(len(a1)):
    lines_a = a1[k].split()
    time3.append(float(lines_a[0]))
    N_lens.append(float(lines_a[1]))

for i in range(len(b1)):
    lines_b = b1[i].split()
    x_ub.append(float(lines_b[0]))
    x_lb.append(float(lines_b[1]))

for j in range(len(c1)):
    lines_c = c1[j].split()
    T.append([float(x) for x in lines_c])

for m in range(len(d1)):
    lines_d = d1[m].split()
    T_ub.append(float(lines_d[0]))
    Tl.append(float(lines_d[1]))

for n in range(len(f1)):
    lines_f = f1[n].split()
    xf3.append(float(lines_f[1]))

for p in range(len(g1)):
    lines_g = g1[p].split()
    time2.append(float(lines_g[0]))
    xf2.append(float(lines_g[1]))

for q in range(len(h1)):
    lines_h = h1[q].split()
    time1.append(float(lines_h[0]))

for r in range(len(l1)):
    lines_l = l1[r].split()
    x_dry.append(float(lines_l[1]))

time3 = np.array(time3)
time2 = np.array(time2)
time1 = np.array(time1)
N_lens = np.array(N_lens)
Tl = np.array(Tl) - 273.15
T_ub = np.array(T_ub) - 273.15
x_ub = np.array(x_ub)
x_lb = np.array(x_lb)
xf2 = np.array(xf2)
xf3 = np.array(xf3)
x_dry = np.array(x_dry)
Tf = -0.054 * 2 / 3
T1 = T[:len(time1)]
T2 = T[len(time1) : len(time1) + len(time2)]
T3 = T[len(time1) + len(time2):]

x = np.linspace(0, 1, num = 601)

# 2d plot
fig = pl.figure()
ax = fig.add_subplot(111)

for k in range(len(time1)):
    ax.plot(x, time1[k] * np.ones(len(x)), linewidth = 2, color = 'k')

for i in range(len(time2)):
    xt = x
    Tt = T2[i]
    for j in range(len(xt) - 1):
        if (Tf > Tt[j]) & (Tf < Tt[j+1]):
            xt = np.concatenate((xt[:j+1], [xf2[i]], xt[j+1:]))
            Tt = np.concatenate((Tt[:j+1], [Tf], Tt[j+1:]))

    time_t = time2[i] * np.ones(len(xt))
    index = np.where(np.isclose(Tt, Tf))[0][0]
    ax.plot(xt, time_t, linewidth = 2, color = 'k')
    ax.plot(xt[:index], time_t[:index], linewidth = 2, color = 'g')

for m in range(len(time3)):
    xt = x
    Tt = T3[m]
    x_lb_t = []
    x_ub_t = []
    idx_x_lb_t = []
    idx_x_ub_t = []
    for q in range(len(xt) - 1):
        if (Tf > Tt[q]) & (Tf < Tt[q+1]):
            xt = np.concatenate((xt[:q+1], [xf3[m]], xt[q+1:]))
            Tt = np.concatenate((Tt[:q+1], [Tf], Tt[q+1:]))

    for i in range(int(N_lens[:m].sum()), int(N_lens[:m+1].sum())):
        x_lb_t = x_lb[int(N_lens[:m].sum()): int(N_lens[:m+1].sum())]
        x_ub_t = x_ub[int(N_lens[:m].sum()): int(N_lens[:m+1].sum())]
        for j in  range(len(xt) - 1):
            if (x_lb[i] > xt[j]) & (x_lb[i] < xt[j+1]):
                xt = np.concatenate((xt[:j+1], [x_lb[i]], xt[j+1:]))
                Tt = np.concatenate((Tt[:j+1], [Tl[i]], Tt[j+1:]))
                break 
        
        for k in range(len(xt) - 1):
            if (x_ub[i] > xt[k]) & (x_ub[i] < xt[k+1]):
                xt = np.concatenate((xt[:k+1], [x_ub[i]], xt[k+1:]))
                Tt = np.concatenate((Tt[:k+1], [T_ub[i]], Tt[k+1:]))
                break
    for s in range(len(xt) - 1):
        if (xt[s] < x_dry[m]) & (xt[s+1] > x_dry[m]):
            T_dry = Tt[s] + (Tt[s+1] - Tt[s]) * (x_dry[m] - xt[s]) / (xt[s+1] - xt[s])
            xt = np.concatenate((xt[:s+1], [x_dry[m]], xt[s+1:]))
            Tt = np.concatenate((Tt[:s+1], [T_dry], Tt[s+1:]))
            break
    time_t = time3[m] * np.ones(len(xt))

    ax.plot(xt, time_t, linewidth = 2, color = 'k')
    ax.plot(xt[:np.where(np.isclose(Tt, Tf))[0][0] + 1], time_t[:np.where(np.isclose(Tt, Tf))[0][0] + 1], linewidth = 2, color = 'g')
    dry_idx = np.where((xt > x_dry[m]) | (xt == x_dry[m]))[0]
    ax.plot(xt[dry_idx], time_t[dry_idx], linewidth = 2, color = 'yellow')

    for n in range(len(x_lb_t)):
        idx_x_lb_t.append(np.where(np.isclose(xt, x_lb_t[n]))[0][0])
        idx_x_ub_t.append(np.where(np.isclose(xt, x_ub_t[n]))[0][0])
    for r in range(len(idx_x_lb_t)):
        ax.plot(xt[idx_x_ub_t[r]:idx_x_lb_t[r] + 1], time_t[idx_x_ub_t[r]:idx_x_lb_t[r] + 1], linewidth = 2, color = 'b')
        
ax.set_xlabel('Depth (m)', fontsize = 50)
ax.set_ylabel('Time (h)', fontsize = 50)
ax.tick_params(axis='both', labelsize = 50)
# ax.set_ylim(-5, 1)




pl.show()
