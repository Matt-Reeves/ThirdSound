# ThirdSound

Solves for standing wave solutions of third sound (surface waves) on thin-film helium. Uses the numerical method developed by A. J. Roberts. 


The equations for the surface $\eta(x,t)$ and potential $\phi(x,y,t)$ are 

$$ \partial_t \eta + \phi_x \eta_x = \phi_y  $$

$$ \partial_t \phi - \alpha \left[\frac{1}{(h + \eta)^3} - \frac{1}{h^3} \right] + \frac{1}{2}( \phi_x^2 + \phi_y^2) =0  $$


<img width="611" alt="Screen Shot 2024-03-07 at 2 44 27 pm" src="https://github.com/Matt-Reeves/ThirdSound/assets/65841999/9a79b4a0-afb2-4a98-9194-c03124d58ca3">


`doglegSolve('16')` 
