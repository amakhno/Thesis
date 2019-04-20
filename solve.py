import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as pylab
import pandas as pd
from SkyNet import *
import gc
import os

def write_solution(filename ,multiply):
    with open(filename) as f:
        content = f.readlines()

    content = [x.strip() for x in content]

    t_array = []
    sig_array = []

    for index in range(0, len(content)):
        splitted = content[index].rsplit(' ')
        t_array.append(float(splitted[0])/1e9)
        val = float(splitted[1]) * multiply
        if (val > 0):
            sig_array.append(np.log(val))
        else:
            sig_array.append(np.log(np.abs(val)))
        

    def build_a(t_arr):
        if isinstance(t_arr, list):
            if len(t_arr) >= 7:
                a = []
                for index in range(len(t_arr)):
                    t = t_arr[index]
                    a.append([1, t**(-1), t**(-1.0/3), t**(1.0/3), t, t**(5.0/3), np.log(t)])
                return a
            else:
                pass
        else:
            pass

    a = build_a(t_array)
    b = sig_array

    x = []; y = []

    # for index in range(0, 24-7):
    #     x.append(t_array[index])
    #     y.append(np.linalg.solve(a[index:index + 7], b[index:index + 7]))

    # def plot_solution(t_array, solution, block):
    #     if (block != True):
    #         block = False
    #     a1_array = []; b1_array = []
    #     for t_index in range(0, len(t_array)):
    #         current_a1 = 0
    #         for a_index in range(0, len(solution)):
    #             current_a1 += solution[a_index] * a[t_index][a_index]
    #         a1_array.append(current_a1)
    #         b1_array.append(b[t_index])
    #     pylab.figure()
    #     line1, = pylab.plot(t_array, a1_array, label="A", linestyle="--")
    #     line2, = pylab.plot(t_array, b1_array, label="B")    
    #     pylab.legend(handles=[line1, line2], loc=1)
    #     pylab.show(block=block)

    x_lstsq_full, _, _, _ =np.linalg.lstsq(a, b, rcond=1)
    # plot_solution(t_array, x_lstsq, False)
    # -6.461089e+02-2.855431e+02 4.604533e+03-3.632646e+03                      
    #  1.326413e+02-5.591785e+00 2.405325e+03     
    # plot_solution(t_array, [-6.461089e2,-2.855431e2,4.604533e3,-3.632646e3,1.326413e2, -5.591785, 2.405325e3], True)
    # x_lstsq_part, _, _, _ =np.linalg.lstsq(a[2:20], b[2:20], rcond=1)
    # plot_solution(t_array, x_lstsq, True)
    # print(x_lstsq_full, file=f)

    def save_to_file(x_lstsq_full, filename, elem1, elem2):
        save_file = open(filename, 'a')
        print >>save_file,"1"
        print >>save_file, "     {:>5}{:>5}                            cleaw     0.00000e+00          ".format(elem1, elem2)
        prettyfied = []
        for index in range(len(x_lstsq_full)):
            prettyfied.append("{0:1.6e}".format(x_lstsq_full[index]))
        print >>save_file, "{:>13}{:>13}{:>13}{:>13}                      ".format(prettyfied[0], prettyfied[1],prettyfied[2],prettyfied[3])
        print >>save_file, "{:>13}{:>13}{:>13}                                   ".format(prettyfied[4], prettyfied[5],prettyfied[6])
        save_file.close()

    def save_to_file_with_p(x_lstsq_full, filename, elem1, elem2):
        save_file = open(filename, 'a')
        print >>save_file,"5"
        print >>save_file, "         p{:>5}{:>5}    p                  cleaw     0.00000e+00          ".format(elem1, elem2)
        prettyfied = []
        for index in range(len(x_lstsq_full)):
            prettyfied.append("{0:1.6e}".format(x_lstsq_full[index]))
        print >>save_file, "{:>13}{:>13}{:>13}{:>13}                      ".format(prettyfied[0], prettyfied[1],prettyfied[2],prettyfied[3])
        print >>save_file, "{:>13}{:>13}{:>13}                                   ".format(prettyfied[4], prettyfied[5],prettyfied[6])
        save_file.close()

    names = os.path.splitext(os.path.basename(filename))[0].rsplit('-')

    print(x_lstsq_full)
    save_to_file_with_p(x_lstsq_full, "mylib.txt", names[0], names[1])

def calculate_xray(with_my, fileName):
  do_heat = False; do_inv = False; do_screen = False
  with open("X-ray_burst/sunet") as f:
    nuclides = [l.strip() for l in f.readlines()]

  nuclib = NuclideLibrary.CreateFromWinv("X-ray_burst/winvne_v2.0.dat", nuclides)

  opts = NetworkOptions()
  opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
  opts.MassDeviationThreshold = 1.0E-10
  opts.IsSelfHeating = do_heat
  opts.EnableScreening = do_screen
  opts.DisableStdoutOutput = True

  helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

  strongReactionLibrary = REACLIBReactionLibrary("X-ray_burst/reaclib",
    ReactionType.Strong, do_inv, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Strong reactions", nuclib, opts, True, True)
  weakReactionLibrary = REACLIBReactionLibrary("X-ray_burst/reaclib",
    ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Weak reactions", nuclib, opts, True, True)
  symmetricFission = REACLIBReactionLibrary("X-ray_burst/nfis",
    ReactionType.Strong, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Symmetric neutron induced fission with 0 neutrons emitted",
    nuclib, opts, True, True)
  spontaneousFission = REACLIBReactionLibrary("X-ray_burst/sfis",
    ReactionType.Strong, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Spontaneous fission", nuclib, opts, True, True)  
  myLib = REACLIBReactionLibrary("mylib.txt",
    ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "CAP", nuclib, opts, True, True)

  if (with_my):
    reactionLibraries = [strongReactionLibrary, weakReactionLibrary,
      symmetricFission, spontaneousFission, myLib]
  else:
    reactionLibraries = [strongReactionLibrary, weakReactionLibrary,
      symmetricFission, spontaneousFission]  

  screen = SkyNetScreening(nuclib)
  net = ReactionNetwork(nuclib, reactionLibraries, helm, screen, opts)

  #lib = net.GetNuclideLibrary()
  #fout = open("sunet", "w")
  #for n in lib.Names():
    #fout.write("%5s\n" % n)
  #fout.close()

  dat = np.loadtxt("X-ray_burst/traj_skynet")
  density_vs_time = PiecewiseLinearFunction(dat[:,0], dat[:,2], True)
  temperature_vs_time = PiecewiseLinearFunction(dat[:,0], dat[:,1], True)

  Temp0 = dat[0,1]
  t0 = dat[0,0] + 1.0e-20

  if (do_heat):
    tfinal = 1.0E+03
  else:
    tfinal = 1.242595E+03

  # initial abundance
  #nseResult = NSE.CalcFromTemperatureAndDensity(Temp0, density_vs_time(t0),
      #Ye0, net.GetNuclideLibrary())
  #np.savetxt("init_Y", nseResult.Y())

  initY = np.loadtxt("X-ray_burst/init_Y")

  if (do_heat):
    output = net.EvolveSelfHeatingWithInitialTemperature(initY, t0, tfinal,
        Temp0, density_vs_time, fileName)
  else:
    output = net.Evolve(initY, t0, tfinal, temperature_vs_time, density_vs_time,
        fileName)
  output.Close()
  net = None
  output = None
  print("Hello")

def calculate_rprocess(with_my, fileName):
  do_heat = False; do_inv = False; do_screen = False
  with open("r-process/sunet") as f:
    nuclides = [l.strip() for l in f.readlines()]

  # make trajectory
  t0 = 1.0e-3
  tfin = 1.0e9

  nuclib = NuclideLibrary.CreateFromWinv("r-process/winvne_v2.0.dat", nuclides)

  opts = NetworkOptions()
  opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
  opts.MassDeviationThreshold = 1.0E-10
  opts.IsSelfHeating = do_heat
  opts.EnableScreening = do_screen
  opts.DisableStdoutOutput = True

  helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")
  if (with_my):
    strongReactionLibrary = REACLIBReactionLibrary("output",
      ReactionType.Strong, do_inv, LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclib, opts, True, True)
    weakReactionLibrary = REACLIBReactionLibrary("output",
      ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclib, opts, True, True)
  else:
    strongReactionLibrary = REACLIBReactionLibrary("r-process/reaclib",
      ReactionType.Strong, do_inv, LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclib, opts, True, True)
    weakReactionLibrary = REACLIBReactionLibrary("r-process/reaclib",
      ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclib, opts, True, True)
  symmetricFission = REACLIBReactionLibrary("r-process/nfis",
    ReactionType.Strong, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Symmetric neutron induced fission with 0 neutrons emitted",
    nuclib, opts, True, True)
  spontaneousFission = REACLIBReactionLibrary("r-process/sfis",
    ReactionType.Strong, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Spontaneous fission", nuclib, opts, True, True)  
  myLib = REACLIBReactionLibrary("mylib.txt",
    ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "CAP", nuclib, opts, True, True)

  if (with_my):
    reactionLibraries = [strongReactionLibrary, weakReactionLibrary,
      symmetricFission, spontaneousFission, myLib]
  else:
    reactionLibraries = [strongReactionLibrary, weakReactionLibrary,
      symmetricFission, spontaneousFission]  

  screen = SkyNetScreening(nuclib)
  net = ReactionNetwork(nuclib, reactionLibraries, helm, screen, opts)

  #lib = net.GetNuclideLibrary()
  #fout = open("sunet", "w")
  #for n in lib.Names():
    #fout.write("%5s\n" % n)
  #fout.close()

  dat = np.loadtxt("r-process/traj_skynet")
  density_vs_time = PiecewiseLinearFunction(dat[:,0], dat[:,2], True)
  temperature_vs_time = PiecewiseLinearFunction(dat[:,0], dat[:,1], True)

  Ye0 = 0.07
  Temp0 = dat[0,1]
  tfinal = 5.0e8

  # initial abundance
  #nseResult = NSE.CalcFromTemperatureAndDensity(Temp0, density_vs_time(t0),
      #Ye0, net.GetNuclideLibrary())
  #np.savetxt("init_Y", nseResult.Y())

  initY = np.loadtxt("r-process/init_Y")

  if (do_heat):
    output = net.EvolveSelfHeatingWithInitialTemperature(initY, t0, tfinal,
        Temp0, density_vs_time, fileName)
  else:
    output = net.Evolve(initY, t0, tfinal, temperature_vs_time, density_vs_time, fileName)

index = 1
multiply = 1

onlyfiles = [f for f in os.listdir('out-sig') if os.path.isfile(os.path.join('out-sig', f))]

os.remove("mylib.txt")
for file_index in range(0, len(onlyfiles)):
  write_solution('out-sig/' + onlyfiles[file_index], multiply)

calculate_rprocess(True, "out-r/WithMy-" + str(index))
calculate_rprocess(False, "out-r/WithoutMy-" + str(index))

# while index > 30:
#   index = index - 1
#   multiply = 10**(index)
#   os.remove("mylib.txt")
#   for file_index in range(0, len(onlyfiles)):
#     write_solution('out-sig/' + onlyfiles[file_index], multiply)
#   calculate_rprocess(True, "out-r/WithMy-" + str(index))
#   calculate_rprocess(False, "out-r/WithoutMy-" + str(index))