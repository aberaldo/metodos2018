function sis = f(t,x);
 sis(1) = x(3);
 sis(2) = x(4);
 sis(3) = (-2*x(2)*x(3)*x(4))/(x(1)^2+x(2)^2+1);
 sis(4) = (-2*x(1)*x(3)*x(4))/(x(1)^2+x(2)^2+1);  
endfunction

function [uOut,vOut,pOut,qOut]=paso_euler_adelante(u,v,p,q,h)	
uOut= p*h+u;
 vOut= q*h+v;
 pOut= (-2*v*p*q*h)/(u^2+v^2+1)+p;
 qOut= (-2*u*p*q*h)/(u^2+v^2+1)+q;
end


function [u_n1_k1,v_n1_k1,p_n1_k1,q_n1_k1]=paso_trapecio(u_n,v_n,p_n,q_n,h,toler)	
	
  k=0;
  u_n1_k= u_n+h*p_n;
  v_n1_k= v_n+h*q_n;
  p_n1_k= p_n+h*(-2*v_n*p_n*q_n/(u_n^2+v_n^2+1));
  q_n1_k= q_n+h*(-2*u_n*p_n*q_n/(u_n^2+v_n^2+1));
  
  u_n1_k1= 0;
  v_n1_k1= 0;
  p_n1_k1= 0;
  q_n1_k1= 0;
   
  while (abs(u_n1_k1-u_n1_k) >= toler && abs(v_n1_k1-v_n1_k) >= toler && abs(p_n1_k1-p_n1_k) >= toler && abs(q_n1_k1-q_n1_k) >= toler && k < 10)
 
    k= k+1;
    u_n1_k1= u_n+(h/2)*(p_n+p_n1_k);
    v_n1_k1= v_n+(h/2)*(q_n+q_n1_k);
    p_n1_k1= p_n+(h/2)*((-2*v_n*p_n*q_n)/(u_n^2+v_n^2+1)+(-2*v_n1_k*p_n1_k*q_n1_k)/(u_n1_k^2+v_n1_k^2+1));
    q_n1_k1= q_n+(h/2)*((-2*u_n*p_n*q_n)/(u_n^2+v_n^2+1)+(-2*u_n1_k*p_n1_k*q_n1_k)/(u_n1_k^2+v_n1_k^2+1));
       
  endwhile

end
function [u_n1_k1,v_n1_k1,p_n1_k1,q_n1_k1]=paso_euler_atras(u_n,v_n,p_n,q_n,h,toler)	
	
  k=0;
  u_n1_k= u_n+h*p_n;
  v_n1_k= v_n+h*q_n;
  p_n1_k= p_n+h*(-2*v_n*p_n*q_n/(u_n^2+v_n^2+1));
  q_n1_k= q_n+h*(-2*u_n*p_n*q_n/(u_n^2+v_n^2+1));
  
  u_n1_k1= 0;
  v_n1_k1= 0;
  p_n1_k1= 0;
  q_n1_k1= 0;
   
  while (abs(u_n1_k1-u_n1_k) >= toler && abs(v_n1_k1-v_n1_k) >= toler && abs(p_n1_k1-p_n1_k) >= toler && abs(q_n1_k1-q_n1_k) >= toler && k < 10)
 
    k= k+1;
    u_n1_k1= u_n+(h)*(p_n1_k);
    v_n1_k1= v_n+(h)*(q_n1_k);
    p_n1_k1= p_n+(h)*(-2*v_n1_k*p_n1_k*q_n1_k)/(u_n1_k^2+v_n1_k^2+1);
    q_n1_k1= q_n+(h)*(-2*u_n1_k*p_n1_k*q_n1_k)/(u_n1_k^2+v_n1_k^2+1);
       
  endwhile

end
function [uSig,vSig,pSig,qSig]=paso_heun(u,v,p,q,h,toler)	
	
  divisor= u^2+v^2+1;
  
  uAux= u+h*p;
  vAux= v+h*q;
  pAux= p+h*(-2*v*p*q/divisor);
  qAux= q+h*(-2*u*p*q/divisor);
  
  divisor2= uAux^2+vAux^2+1;
  
  uSig= u+(h/2)*(p+pAux);
  vSig= v+(h/2)*(q+qAux);
  pSig= p+(h/2)*((-2*v*p*q/divisor)+(-2*vAux*pAux*qAux/divisor2));
  qSig= q+(h/2)*((-2*u*p*q/divisor)+(-2*uAux*pAux*qAux/divisor2));
end

function [uSig,vSig,pSig,qSig]=paso_RK(u,v,p,q,h)	

  divisor= u^2+v^2+1;

	k1_1 = p;
	k1_2 = q;
  k1_3 = -2*v*p*q/divisor;
  k1_4 = -2*u*p*q/divisor;
  
  uAux = u+(1/2)*k1_1*h;
  vAux = v+(1/2)*k1_2*h;
  pAux = p+(1/2)*k1_3*h;
  qAux = q+(1/2)*k1_4*h;
  
  divisor= uAux^2+vAux^2+1;

	k2_1 = pAux;
	k2_2 = qAux;
  k2_3 = -2*vAux*pAux*qAux/divisor;
  k2_4 = -2*uAux*pAux*qAux/divisor;
  
  uAux = u+(1/2)*k2_1*h;
  vAux = v+(1/2)*k2_2*h;
  pAux = p+(1/2)*k2_3*h;
  qAux = q+(1/2)*k2_4*h;
  
  divisor= uAux^2+vAux^2+1;

  k3_1 = pAux;
	k3_2 = qAux;
  k3_3 = -2*vAux*pAux*qAux/divisor;
  k3_4 = -2*uAux*pAux*qAux/divisor;
  
  uAux = u+k3_1*h;
  vAux = v+k3_2*h;
  pAux = p+k3_3*h;
  qAux = q+k3_4*h;

  divisor= uAux^2+vAux^2+1;

  k4_1 = pAux;
	k4_2 = qAux;
  k4_3 = -2*vAux*pAux*qAux/divisor;
  k4_4 = -2*uAux*pAux*qAux/divisor;
  
  
	uSig= u+(h/6)*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1);

	vSig= v+(h/6)*(k1_2 + 2*k2_2  + 2*k3_2  + k4_2 );
	
	pSig = p +(h/6)*(k1_3 + 2*k2_3 + 2*k3_3 + k4_3);
	
	qSig = q +(h/6)*(k1_4 + 2*k2_4 + 2*k3_4 + k4_4);
	
end

t = (0:0.01:1);

#ode45#
[odeT,odeRes] = ode45("f", t', [0,0,1,2]);

#EULER HACIA ADELANTE#
toler= 10^-4;
h= 0.01;
u= 0; v= 0; p= 1; q= 2;
RES_EA(1,:)=[u,v,p,q];

for i=0:100
  [u,v,p,q]=paso_euler_adelante(u,v,p,q,h);
  RES_EA(i+1,:)=[u,v,p,q];
end

#TRAPECIO#
toler= 10^-4;
h= 0.01;
u= 0; v= 0; p= 1; q= 2;
RES_T(1,:)=[u,v,p,q];

for i=0:100
  [u,v,p,q]=paso_trapecio(u,v,p,q,h,toler);
  RES_T(i+1,:)=[u,v,p,q];
end  

#EULER HACIA ATRAS#
toler= 10^-4;
h= 0.01;
u= 0; v= 0; p= 1; q= 2;
RES_AT(1,:)=[u,v,p,q];

for i=0:100
  [u,v,p,q]=paso_euler_atras(u,v,p,q,h,toler);
  RES_AT(i+1,:)=[u,v,p,q];
end

#HEUN#
toler= 10^-4;
h= 0.01;
u= 0; v= 0; p= 1; q= 2;
RES_H(1,:)=[u,v,p,q];

for i=0:100
  [u,v,p,q]=paso_heun(u,v,p,q,h,toler);
  RES_H(i+1,:)=[u,v,p,q];
end

#RK#
toler= 10^-4;
h= 0.01;
u= 0; v= 0; p= 1; q= 2;
RES_RK(1,:)=[u,v,p,q];

for i=0:100
  [u,v,p,q]=paso_RK(u,v,p,q,h,toler);
  RES_RK(i+1,:)=[u,v,p,q];
end

#ADAMS BASHFORTH#

function [uOut,vOut,pOut,qOut]=paso_AB_2(u_1, v_1, p_1, q_1, u, v, p, q, h)	
	uOut= u_1 + h*((3/2)*p_1 - (1/2)*p);
	vOut= v_1 + h*((3/2)*q_1 - (1/2)*q);
	pOut= p_1 + h*((3/2)*((-2*v_1*p_1*q_1)/(u_1^2+v_1^2+1)) - (1/2)*(-2*v*p*q)/(u^2+v^2+1));
	qOut= q_1 + h*((3/2)*((-2*u_1*p_1*q_1)/(u_1^2+v_1^2+1)) - (1/2)*(-2*u*p*q)/(u^2+v^2+1));
end

toler= 10^-4;
h= 0.01;
u= 0; v= 0; p= 1; q= 2;
[u_1,v_1,p_1,q_1]=paso_euler_adelante(u,v,p,q,h);

AB(1,:)=[u,v,p,q];
AB(2,:)=[u_1,v_1,p_1,q_1];

for i=2:100
	
	[u_2,v_2,p_2,q_2]=paso_AB_2(u_1,v_1,p_1,q_1,u,v,p,q,h);
	AB(i+1,:)=[u_2,v_2,p_2,q_2];
	
	u = u_1;
	v = v_1;
	p = p_1;
	q = q_1;
	
	u_1 = u_2;
	v_1 = v_2;
	p_1 = p_2;
	q_1 = q_2;
end

#GRAFICA u#
figure(1);
plot (t,RES_EA(:,1),";Euler hacia adelante;", t,RES_T(:,1),";Trapecio;", t,RES_AT(:,1),";Euler hacia atras;", t,RES_H(:,1),";Heun;", t,RES_RK(:,1),";Runge-Kutta;", t,odeRes(:,1),";Ode45;", t,AB(:,1),";AB;");
title ("Grafica u");
xlabel ("h","fontsize", 14);
ylabel ("f","fontsize", 14);

#GRAFICA v#
figure(2);
plot (t,RES_EA(:,2),";Euler hacia adelante;", t,RES_T(:,2),";Trapecio;", t,RES_AT(:,2),";Euler hacia atras;", t,RES_H(:,2),";Heun;", t,RES_RK(:,2),";Runge-Kutta;", t,odeRes(:,2),";Ode45;", t,AB(:,2),";AB;");
title ("Grafica v");
xlabel ("h","fontsize", 14);
ylabel ("f","fontsize", 14);

#GRAFICA p#
figure(3);
plot (t,RES_EA(:,3),";Euler hacia adelante;", t,RES_T(:,3),";Trapecio;", t,RES_AT(:,3),";Euler hacia atras;", t,RES_H(:,3),";Heun;", t,RES_RK(:,3),";Runge-Kutta;", t,odeRes(:,3),";Ode45;", t,AB(:,3),";AB;");
title ("Grafica p");
xlabel ("h","fontsize", 14);
ylabel ("f","fontsize", 14);

#GRAFICA q#
figure(4);
plot (t,RES_EA(:,4),";Euler hacia adelante;", t,RES_T(:,4),";Trapecio;", t,RES_AT(:,4),";Euler hacia atras;", t,RES_H(:,4),";Heun;", t,RES_RK(:,4),";Runge-Kutta;", t,odeRes(:,4),";Ode45;", t,AB(:,4),";AB;");
title ("Grafica q");
xlabel ("h","fontsize", 14);
ylabel ("f","fontsize", 14);