#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

template<typename T>
class UndirectedGraph {
private:
	std::map<T, std::vector<T>> list;
	void AddPair(T u, T v) {
		auto i = list.find(u);
		if (i == list.end()) {
			std::vector<T> temp;
			temp.push_back(v);
			list[u] = temp;
		}
		else {
			(i->second).push_back(v);
		}
	}
public:
	void AddEdge(T u, T v) {
		AddPair(u, v);
		AddPair(v, u);
	}
	void Vertex(T u) {
		auto i = list.find(u);
		/*while (i < list.end() && *(*i).begin() != u) {
			i++;
		}*/
		if (i != list.end()) {
			for (auto j = (i->second).begin(); j < (i->second).end(); j++) {
				std::cout << *j << " ";
			}
		}
		std::cout << "\n";
	}
};

int main() {
	UndirectedGraph<int> graph;
	int CommandsNumber, VertexNumber;
	std::cin >> VertexNumber >> CommandsNumber;

	for (int i = 0; i < CommandsNumber; i++) {
		int command;
		std::cin >> command;
		if (command == 1) {
			int u, v;
			std::cin >> u >> v;
			graph.AddEdge(u, v);
		}
		if (command == 2) {
			int u;
			std::cin >> u;
			graph.Vertex(u);
		}
	}
	//system("pause");
	return 0;
}