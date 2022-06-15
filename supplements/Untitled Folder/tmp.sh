lstH = []

with open('compHawk') as f:
	for line in f:
		lstH.append(line.strip())

with open('comp') as f:
	for line in f:
		line = line.strip()
		print(line)
		if line in lstH:
			print(line)
		


