import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import os

AINDEX = 97

class Graph:
    # Конструктор класса
    def __init__(self, size, bFullGraph = False, bMultiedge = False, bLoop = False):
        self._nodes = size
        self._bFullGraph = bFullGraph
        self._bMultiedge = bMultiedge
        self._bLoop = bLoop
        self._imatrix = None
        self._amatrix = list()
        for i in range(self._nodes):
            self._amatrix.append(list())
            for j in range(self._nodes):
                self._amatrix[i].append(0)
    
    # Вывод матрицы смежности
    def showAdjacencyMatrix(self):
        print("Матрица смежности: ")
        for i in range(self._nodes * 2 + 1):
            for j in range(self._nodes * 2 + 1):
                if j % 2 == 1:
                    print(end = "  ") # |
                elif i % 2 == 1:
                    print(end = " ") # -
                elif (i == 0 and j != 0) or (j == 0 and i != 0):
                    print(end = f"{chr(AINDEX + (i + j) // 2 - 1)}")
                elif i // 2 > 0 and j // 2 > 0:
                    print(end = f"{self._amatrix[i // 2 - 1][j // 2 - 1]}")
                else:
                    print(end = " ")
            print()

    # Вывод графа (при помощи networkx)
    def showGraph(self):
        nodeMap = dict()
        edgeMap = dict()
        loopMap = dict()
        for i in range(0, self._nodes):
            nodeMap.update({i: chr(AINDEX + i)})

        count = 1
        for i in range(self._nodes):
            for j in range(i, self._nodes):
                if self._amatrix[i][j] != 0:
                    edgeName = ""
                    if self._amatrix[i][j] > 1:
                        for k in range(self._amatrix[i][j] - 1):
                            edgeName += f"{count}, "
                            count += 1
                        edgeName += f"{count}"
                        count += 1
                    else:
                        edgeName = f"{count}"
                        count += 1
                    if i == j:
                        edgeName += "\n\n"
                        loopMap.update({(chr(AINDEX + i), chr(AINDEX + j)): edgeName})
                    else:
                        edgeMap.update({(chr(AINDEX + i), chr(AINDEX + j)): edgeName})

        G = nx.Graph(np.array(self._amatrix))
        nx.relabel_nodes(G, nodeMap, False)
        pos = nx.circular_layout(G)
        nx.draw(G, pos, with_labels = True)
        nx.draw_networkx_edge_labels(G, pos, edge_labels = edgeMap)
        nx.draw_networkx_edge_labels(G, pos, edge_labels = loopMap)
        plt.show()

    # Заполнение матрицы смежности при помощи рандома
    def setRandomMatrix(self):
        deltaIndex = 0 if self._bLoop else 1
        minEdges = 1 if self._bFullGraph else 0
        maxEdges = 3 if self._bMultiedge else 1

        for i in range(0, self._nodes):
            for j in range(i + deltaIndex, self._nodes):
                self._amatrix[i][j] = self._amatrix[j][i] = random.randint(minEdges, maxEdges) 
    
    # Подсчет ребер графа
    def updateEdges(self):
        edges = 0
        for i in range(self._nodes):
            for j in range(i, self._nodes):
                edges += self._amatrix[i][j]
        self._edges = edges
    
    # Создание матрицы инцидентности
    def makeIncidenceMatrix(self):
        self.updateEdges()

        self._imatrix = list()
        for i in range(self._nodes):
            self._imatrix.append(list())
            for j in range(self._edges):
                self._imatrix[i].append(0)
        
        edgeIndex = 0
        for i in range(self._nodes):
            for j in range(i, self._nodes):
                if self._amatrix[i][j] != 0:
                    value = 1
                    if i == j:
                        value = 2
                    
                    for k in range(self._amatrix[i][j]):
                            self._imatrix[i][edgeIndex] = self._imatrix[j][edgeIndex] = value
                            edgeIndex += 1

    # Вывод матрицы инцидентности
    def showIncidenceMatrix(self):
        if self._imatrix is None:
            self.makeIncidenceMatrix()
        print("Матрица инцидентности: ")
        for i in range(self._nodes * 2 + 1):
            for j in range(self._edges * 2 + 1):
                if j % 2 == 1:
                    print(end = "  ") # |
                elif i % 2 == 1:
                    print(end = " ") # -
                elif j == 0 and i != 0:
                    print(end = f"{chr(AINDEX + (i + j) // 2 - 1)}")
                elif i == 0 and j != 0:
                    print(end = f"{j // 2}")
                elif i // 2 > 0 and j // 2 > 0:
                    print(end = f"{self._imatrix[i // 2 - 1][j // 2 - 1]}")
                    if j >= 20:
                        print(end = " ")
                else:
                    print(end = " ")
            print()
                    



def main():
    menu = "Выберите режим работы программы (можно выбрать несколько):\n"
    menu += "Enter - Обычный режим\n1 - Полный граф\n2 - Кратные ребра\n3 - Наличие петель (могут быть кратными)\nq - Выход\n"

    mode = input(menu)

    while ('q' not in mode):
        NodeNumber = int(input("Введите количество вершин графа: "))
        graph = Graph(NodeNumber, '1' in mode, '2' in mode, '3' in mode)
        graph.setRandomMatrix()
        graph.showAdjacencyMatrix()
        print()
        graph.showIncidenceMatrix()
        graph.showGraph()
        os.system("cls")
        mode = input(menu)
    os.system("cls")


    


if __name__ == "__main__":
    main()