from random import randint
import subprocess

for n in range(0, 10):

	particles = "100 1e4 "
	atomic = "./sol6atom "
	reduction = "./sol6reduc "
	serial = "./sol2serial "
	pos = ""

	for i in range(n):
		for j in range(n):
			for k in range(n):
				pos = str(i * 5e-9) + " " + str(j * 5e-9) + " " + str(k * 5e-9) + " "
				velx = vely = velz = "0 "

				# With probability 1/5, give the components a starting velocity
				if(randint(1,10) > 8): 
					velx = "1e-12 "

				if(randint(1,10) > 8):
					vely = "1e-12 "

				if(randint(1,10) > 8):
					vely = "1e-12 "

				particles += pos + velx + vely + velz + "39.948 "

	atomic += particles
	reduction += particles
	serial += particles

	with open('AtomicTimings.txt', "a") as outfile:
		subprocess.run(atomic.split(), stdout=outfile)

	with open('ReductionTimings.txt', "a") as outfile:
		subprocess.run(reduction.split(), stdout=outfile)

	with open('SerialTimings.txt', "a") as outfile:
		subprocess.run(serial.split(), stdout=outfile)
