SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 6, 6, 6};
//+
Physical Surface("sup1", 13) = {6};
//+
Physical Surface("sup2", 14) = {1};
//+
Physical Surface("sup3", 15) = {5};
//+
Physical Surface("sup4", 16) = {2};
//+
Physical Surface("sup6", 17) = {4};
//+
Physical Surface("sup5", 18) = {3};
//+
Physical Volume("vol", 19) = {1};

// --- CONFIGURAÇÕES DA MALHA ---
Mesh.RecombineAll = 1;       // Transforma triângulos em quadriláteros
Mesh.Algorithm = 8;          // Algoritmo "Frontal-Delaunay para quadriláteros"
Mesh.ElementOrder = 2;       // Define segunda ordem (quadrática)
Mesh.SecondOrderLinear = 0;  // Ativa os nós centrais (Lagrangeano / 9 nós)
