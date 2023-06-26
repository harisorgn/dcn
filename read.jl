function read_genesis(dir_path)
    file = open("./data/pDend.txt", "r")
	pdend = readdlm(file) ;
	close(file)

	file = open("./data/dDend.txt", "r")
	ddend = readdlm(file) ;
	close(file)

	file = open("./data/axIS.txt", "r")
	axis = readdlm(file) ;
	close(file)

	file = open("./data/axIN.txt", "r")
	axin = readdlm(file) ;
	close(file)

    return (pdend, ddend, axis, axin)
end