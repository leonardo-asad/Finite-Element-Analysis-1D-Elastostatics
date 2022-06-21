//Standard C++ libraries
#include <iostream>

//Find the value of xi at the given node (using deal.II node numbering)
int basisFunctionOrder = 3;

double xi_at_node(unsigned int dealNode){
  double xi;

  if(dealNode == 0){
    xi = -1.;
  }
  else if(dealNode == 1){
    xi = 1.;
  }
  else if(dealNode <= basisFunctionOrder){
    xi = -1. + 2.*(dealNode-1.)/basisFunctionOrder;
  }
  else{
    std::cout << "Error: you input node number "
	      << dealNode << " but there are only "
	      << basisFunctionOrder + 1 << " nodes in an element.\n";
    exit(0);
  }

  return xi;
}

//Define basis functions
double basis_function(unsigned int node, double xi){
  /*"basisFunctionOrder" defines the polynomial order of the basis function,
    "node" specifies which node the basis function corresponds to,
    "xi" is the point (in the bi-unit domain) where the function is being evaluated.
    You need to calculate the value of the specified basis function and order at the given quadrature pt.*/

  double value = 1.; //Store the value of the basis function in this variable

  /*You can use the function "xi_at_node" (defined above) to get the value of xi (in the bi-unit domain)
    at any node in the element - using deal.II's element node numbering pattern.*/

  for (int i = 0; i < basisFunctionOrder + 1; i++){
    if (i != node){
      value *= (xi - xi_at_node(i))/(xi_at_node(node) - xi_at_node(i));
    }
  }

  return value;
}

//Define basis function gradient
double basis_gradient(unsigned int node, double xi){
  /*"basisFunctionOrder" defines the polynomial order of the basis function,
    "node" specifies which node the basis function corresponds to,
    "xi" is the point (in the bi-unit domain) where the function is being evaluated.
    You need to calculate the value of the derivative of the specified basis function and order at the given quadrature pt.
    Note that this is the derivative with respect to xi (not x)*/

  //Store the value of the gradient of the basis function in this variable

  double value = 0.;

  /*You can use the function "xi_at_node" (defined above) to get the value of xi (in the bi-unit domain)
    at any node in the element - using deal.II's element node numbering pattern.*/

  for (int i = 0; i < basisFunctionOrder + 1; i++){
    if (i != node){
      double mul = 1.;
      for (int m = 0; m < basisFunctionOrder + 1; m++){
        if (m != i && m != node){
          mul *= (xi - xi_at_node(m))/(xi_at_node(node) - xi_at_node(m));
        }
      }
      value += (1 / (xi_at_node(node) - xi_at_node(i)))*mul;
    }

  }

  return value;
}



int main() {
  //std::cout << "Hello world";

  //std::cout << xi_at_node(1);

  //std::cout << basis_function(0, 1.0);

  std::cout << basis_gradient(0, 1.0);

  return 0;
}
