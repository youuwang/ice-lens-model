import numpy as np
import matplotlib.pyplot as pl

a = open('t_N.dat', 'r')
b = open('ub_lb_t', 'r')
c = open('results/T_profile.dat', 'r')
d = open('Tu_Tl.dat', 'r')
f = open('results/label_time.dat', 'r')
g = open('t_xf2.dat', 'r')
h = open('t1.dat', 'r')
l = open('t_dry.dat', 'r')
u = open('results/x_time.dat', 'r')

a1 = a.readlines()
b1 = b.readlines()
c1 = c.readlines()
d1 = d.readlines()
f1 = f.readlines()
g1 = g.readlines()
h1 = h.readlines()
l1 = l.readlines()
u1 = u.readlines()

del a
del b
del c
del d
del f
del g
del h
del l
del u

time3 = []
time2 = []
time1 = []
xf2 = []
N_lens = []
Tl = []
T_ub = []
x_ub = []
x_lb = []
T = []
label_t = []
x_dry = []
x_time = []

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
    lines_f = [float(y) for y in lines_f]
    label_t.append(lines_f)

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

for s in range(len(u1)):
    lines_u = u1[s].split()
    lines_u = [float(y) for y in lines_u]
    x_time.append(lines_u)

time3 = np.array(time3)
time2 = np.array(time2)
time1 = np.array(time1)
N_lens = np.array(N_lens)
Tl = np.array(Tl) - 273.15
T_ub = np.array(T_ub) - 273.15
x_ub = np.array(x_ub)
x_lb = np.array(x_lb)
x_dry = np.array(x_dry)
Tf = -0.054
T1 = T[:len(time1)]
T2 = T[len(time1) : len(time1) + len(time2)]
T3 = T[len(time1) + len(time2):]

x = np.linspace(0, 1, num = 601)

# 2d plot
fig = pl.figure()
ax = fig.add_subplot(111)
width = 8

for k in range(len(time1)):
    ax.plot(time1[k] * np.ones(len(x)), x, linewidth = width, color = '0.5')

if len(time2) > 0:
    for i in range(len(time2)):
        xt = x
        Tt = T2[i]
        for j in range(len(xt) - 1):
            if (Tf - Tt[j]) * (Tf - Tt[j+1]) < 0:
                xt = np.concatenate((xt[:j+1], [xf2[i]], xt[j+1:]))
                Tt = np.concatenate((Tt[:j+1], [Tf], Tt[j+1:]))

        time_t = time2[i] * np.ones(len(xt))
        index = np.where(np.isclose(Tt, Tf))[0][0]
        ax.plot(time_t, xt, linewidth = width, color = '0.5')
        ax.plot(time_t[:index], xt[:index], linewidth = width, color = 'g')

for m in range(len(time3)):
    xt = x_time[m]
    Tt = T3[m]
    label_time = label_t[m]


    for q in range(len(xt) - 1):
        if (Tf - Tt[q]) * (Tf - Tt[q+1]) < 0:
            xf = (Tf - Tt[q]) / (Tt[q+1] - Tt[q]) * (xt[q+1] - xt[q]) + xt[q]
            xt = np.concatenate((xt[:q+1], [xf], xt[q+1:]))
            Tt = np.concatenate((Tt[:q+1], [Tf], Tt[q+1:]))
            label_time = np.concatenate((label_time[:q+1], [11], label_time[q+1:]))
        elif Tt[q] == Tf:
            label_time[q] = 11
    

    for r in range(len(xt) - 1):
        if Tt[r] * Tt[r + 1] < 0:
            x0 = Tt[r] / (Tt[r] - Tt[r+1]) * (xt[r+1] - xt[r]) + xt[r]
            xt = np.concatenate((xt[:r+1], [x0], xt[r+1:]))
            Tt = np.concatenate((Tt[:r+1], [0], Tt[r+1:]))
            label_time = np.concatenate((label_time[:r+1], [20], label_time[r+1:]))
        elif Tt[r] == 0:
            label_time[r] = 20

    for i in range(int(N_lens[:m].sum()), int(N_lens[:m+1].sum())):
        for j in  range(len(xt) - 1):
            if (x_lb[i] > xt[j]) & (x_lb[i] < xt[j+1]):
                xt = np.concatenate((xt[:j+1], [x_lb[i]], xt[j+1:]))
                Tt = np.concatenate((Tt[:j+1], [Tl[i]], Tt[j+1:]))
                if Tl[i] > 0:
                    label_time = np.concatenate((label_time[:j+1], [32], label_time[j+1:]))
                
                else:
                    label_time = np.concatenate((label_time[:j+1], [22], label_time[j+1:]))
            elif x_lb[i] == xt[j]:
                if (Tl[i] < 0) | (Tl[i] == 0):
                    label_time[j] = 22
                else:
                    label_time[j] = 32

            
        for k in range(len(xt) - 1):
            if (x_ub[i] > xt[k]) & (x_ub[i] < xt[k+1]):
                xt = np.concatenate((xt[:k+1], [x_ub[i]], xt[k+1:]))
                Tt = np.concatenate((Tt[:k+1], [T_ub[i]], Tt[k+1:]))
                if T_ub[i] > 0:
                    label_time = np.concatenate((label_time[:k+1], [31], label_time[k+1:]))
                
                else:
                    label_time = np.concatenate((label_time[:k+1], [21], label_time[k+1:]))
            elif x_ub[i] == xt[k]:
                if (T_ub[i] < 0) | (T_ub[i] == 0):
                    label_time[k] = 21
                else:
                    label_time[k] = 31

    for s in range(len(xt) - 1):
        if (xt[s] < x_dry[m]) & (xt[s+1] > x_dry[m]):
            T_dry = Tt[s] + (Tt[s+1] - Tt[s]) * (x_dry[m] - xt[s]) / (xt[s+1] - xt[s])
            xt = np.concatenate((xt[:s+1], [x_dry[m]], xt[s+1:]))
            Tt = np.concatenate((Tt[:s+1], [T_dry], Tt[s+1:]))
            label_time = np.concatenate((label_time[:s+1], [1], label_time[s+1:]))
            break
    
    time_t = time3[m] * np.ones(len(xt))

    for n in range(len(xt) - 1):
        if (label_time[n] == 1) & (Tt[n] < Tf):
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'g')
        elif label_time[n] == 2:
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'b')
        elif label_time[n] == 3:
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'r')
        elif (label_time[n] == 11) & (Tt[n+1] < Tf) & ((label_time[n+1] == 1) | (label_time[n+1] == 21) | (label_time[n+1] == 31)):
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'g')
        elif (label_time[n] == 11) & (label_time[n+1] == 2):
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'b')
        elif (label_time[n] == 11) & (label_time[n+1] == 20):
            if n < len(xt) - 2:
                if (label_time[n+2] == 3) | ((label_time[n+2] == 32) | (label_time[n+2] == 22)):
                    ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'b')
                else:
                    ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = '0.5')
        elif (label_time[n] == 11) & (label_time[n+1] == 22):
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'b')
        elif (label_time[n] == 20) & ((label_time[n+1] == 2) | (label_time[n+1] == 22)):
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'b')
        elif (label_time[n] == 20) & (label_time[n+1] == 11):
            if n < len(xt) - 2:
                if (label_time[n+2] == 2) | (label_time[n+2] == 22):
                    ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'b')
                else:
                    ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = '0.5')
        elif (label_time[n] == 20) & ((label_time[n+1] == 3) | (label_time[n+1] == 32)):
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'r')
        elif label_time[n] == 21:
            if (Tt[n + 1] > 0):
                ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'r')
            else:
                ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'b')
        elif label_time[n] == 22:
            if Tt[n] < Tf:
                ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'g')
            else:
                ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = '0.5')
        elif label_time[n] == 31:
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = 'r')
        else:
            ax.plot(time_t[n:n+2], xt[n:n+2], linewidth = width, color = '0.5')
    
    dry_idx = np.where((xt > x_dry[m]) | (xt == x_dry[m]))[0]
    ax.plot(time_t[dry_idx], xt[dry_idx], linewidth = width, color = 'yellow')

ax.set_ylabel('Depth (m)', fontsize = 60, labelpad = 20)
ax.set_xlabel('Time (h)', fontsize = 60, labelpad = 20)
ax.tick_params(axis='both', labelsize=60, pad=20)
# ax.set_ylim(0, 0.4)
# ax.set_xlim(30, 20)
# ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])

pl.gca().set_aspect(10/0.4, adjustable='box')
pl.subplots_adjust(top = 0.95, bottom = 0.2)

pl.show()
