##########				     ##############
########## SKRYPT DZIALA TYLKO NA PYTHON 2.7 ##############
##########   KONIECZNE JEST ZAINSTALOWANIE   ##############
##########      CABSflex STANDALONE	     ##############
##########                                   ##############

dir_path = 'project12520' #### Nalezy podmienic 'pdb_files' na sciezke do katalogu z plikami pdb,####
			  ####   ktore ma analizowac CABSflex, a nastepnie skompilowac skrypt.   ####

		#### Przy kompilacji, program wczytuje pliki z katalogu dir_path i kolejno puszcza 
		#### CABSa na tych plikach. Nastepnie sprawdza energie poszczegolnych konformacji
		#### i w katalogu CABS_wyniki (w pliku best_10) znajduja sie sciezki do 10 
		#### (roznych - tzn. wzietych z roznych wejsciowych pdb) modeli o najmniejszych energiach.


import subprocess
import multiprocessing
import os
import csv


def CABSflex(dir_path, file_path, counter, q):

	subprocess.call('CABSflex -i ' + dir_path + '/' + file_path + ' --work-dir ./CABS_wyniki/' + str(counter) + file_path, shell = True) #korzystam z Windowsa, ale wydaje mi sie ze pod linuxem tez powinno zadzialac
	
	plot_paths = os.listdir('CABS_wyniki/' + str(counter) + file_path + '/plots')
	#print(plot_paths)
	
	csv_path = plot_paths[0]#plot_paths[0] to plik csv z rmsd i enegia kazdej konformacji

	f = open('CABS_wyniki/' + str(counter) + file_path + '/plots/' + csv_path, 'r') #otwieram ten plik csv i tworze slownik gdzie rmsd jest kluczem a odpowiadajace mu energie sa wartosciami
	reader = csv.reader(f)
	models = {}
	key_set = set()
	for row in reader:
		data = row[0].split('\t')
		i = 0
		key = data[0] + '_(' + str(i) + ')'
		while key in key_set:
			i += 1
			key = data[0] + '_(' + str(i) + ')'
		models[key] = data[1]
		key_set.add(key)
	f.close()


	output_data_paths = os.listdir('CABS_wyniki/' + str(counter) + file_path + '/output_data') 
	txt_path = output_data_paths[2] #w output_data_paths[2] jest plik txt (medoids_rmsds) z rmsd dla 10 przedstawicieli klastrow. Bede uzywal tych rmsd zeby znalezc energie odpowiadajace danym modelom 
	f2 = open('CABS_wyniki/' + str(counter) + file_path + '/output_data/' + txt_path, 'r')

	E_medoids = {} #klauczami sa rmsd z pliku medoids_rmsds, a wartosciami energie odpowiadajace tym rmsd z pliku csv
	for i in range(10):
		E_medoids[i] = [] #wartosciami slownika E_medoids sa listy, bo po zaogrogleniu rmsd z medoids_rmsds do 3 miejsc po przecinku co jest konieczne do porownywania ich z rmsd z pliku scv (slownik models), jedno rmsd pasuje do kilku energi

	nr_modelu = 0
	for rmsd in f2.read()[:-2].strip().split(';\n'):
		i = 0
		#print(rmsd)
		#print('teraz ostatnie rmsd')
		key = str(round(float(rmsd),3)) + '_(' + str(i) + ')'
		#print(key)
		while key in key_set:
			E_medoids[nr_modelu].append(float(models[key]))
			i += 1
			key = str(round(float(rmsd),3)) + '_(' + str(i) + ')'
			#print(key)
		nr_modelu += 1
	print(E_medoids)
	f2.close()

	minimalna = 0

	E_best_medoids = {} #W tym slowniku bede przechowywal energie jako klucz i sciezke do moedelu jako wartosc
	print(E_medoids)
	for i in range(10): #dla kazdego rmsd z E_medoids wybieram najlepsza energie
		if len(E_medoids[i]) > 0:
			mini = min(E_medoids[i])
		else:
			print("pusty klucz rmsd")
			mini = 0
		if mini < minimalna:
			model_path = 'CABS_wyniki/' + str(counter) + file_path + '/output_pdbs' + 'model_' + str(i)
			minimalna = mini
		E_best_medoids[i] = mini
	#print(E_best_medoids)

	with open('CABS_wyniki/' + str(counter) + file_path + '/best_E.csv', 'w') as f3:
	    for key in E_best_medoids.keys():
	        f3.write("%s,%s\n"%(key,E_best_medoids[key]))
	
	q.put([minimalna, model_path])

pdb_paths = os.listdir(dir_path)
#print (pdb_paths)
    
def QuickSort(A, l=0, r=None):
	if r - l > 1:
		q = Partition(A,l,r)
		QuickSort(A,l,q)
		QuickSort(A,q+1,r)

def Partition(A,l,r):
	i = l + 1
	j = r - 1
	while True:
		while i < r and A[i][0] <= A[l][0]:
			i += 1
		while A[j][0] > A[l][0]:
			j -= 1
		if i < j:
			(A[i], A[j]) = (A[j],A[i])
		else:
			break
	(A[l], A[j]) = (A[j],A[l])
	return j


if __name__ == '__main__':
	
	processes = []
	q = multiprocessing.Queue()
	best_models = []
	
	for counter, file_path in enumerate(pdb_paths): #biore po kolei kazdy plik z dir_path, przepuszczam przez CABSa i wyniki zapisuje w katalogu CABS_wyniki
		p = multiprocessing.Process(target = CABSflex, args = [dir_path, file_path, counter, q])
		p.start()
		processes.append(p)
	
	for process in processes:
		process.join()

	while not q.empty():
		best_models.append(q.get())

	QuickSort(best_models, r = len(best_models))

	with open('CABS_wyniki/best_10_E.csv', 'w') as f4:
		for model in best_models:
		    f4.write("%s,%s\n"%(model[0],model[1]))





















