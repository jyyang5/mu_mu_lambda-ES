# mu_mu_lambda-ES

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

Implement a (mu/mu,lambda)-ES

The algorithm is initialized with mu parents, choose the centroid of parents to generate lambda offsprings and then pick the mu best offspring as the parents of next iteration.

**Note**: lambda >= mu

## Use CSA to adapt step size

Replace the previous step size adaptation using sigma = sigma_star/n*dist with Cumulative Step Length Adaptation (CSA) from https://pdfs.semanticscholar.org/2556/e4670234fc74232a1a419260fe8c1a29d0ba.pdf.
$c = frac{}{sqrt{n}}$