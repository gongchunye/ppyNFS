def loadParamsFile():
	params = {}
	with open("params.txt") as paramsFile:
		for line in paramsFile:
			name, var = line.partition("=")[::2]
			params[name.strip()] = var
	return (int(params["n"]),eval(params["nfsPoly"]),int(params["m"]),int(params["B"]),int(params["M"]),int(params["K"]))

def loadFileArray(fname):
	array = []
	with open(fname) as file:
		for line in file:
			array.append(eval(line))
			
	return array