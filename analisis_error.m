
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
   
  while (abs(u_n1_k1-u_n1_k) >= toler && abs(v_n1_k1-v_n1_k) >= toler && abs(p_n1_k1-p_n1_k) >= toler && abs(q_n1_k1-q_n1_k) >= toler  && k<10)
 
    k= k+1;
    u_n1_k1= u_n+(h/2)*(p_n+p_n1_k);
    v_n1_k1= v_n+(h/2)*(q_n+q_n1_k);
    p_n1_k1= p_n+(h/2)*((-2*v_n*p_n*q_n)/(u_n^2+v_n^2+1)+(-2*v_n1_k*p_n1_k*q_n1_k)/(u_n1_k^2+v_n1_k^2+1));
    q_n1_k1= q_n+(h/2)*((-2*u_n*p_n*q_n)/(u_n^2+v_n^2+1)+(-2*u_n1_k*p_n1_k*q_n1_k)/(u_n1_k^2+v_n1_k^2+1));
       
  endwhile

end
function [u_n1_k1,v_n1_k1,p_n1_k1,q_n1_k1]=paso_euler_atras(u_n,v_n,p_n,q_n,h,toler)	
	
  #fprintf('euler atras paso \n');
  k=0;
  u_n1_k= u_n+h*p_n;
  v_n1_k= v_n+h*q_n;
  p_n1_k= p_n+h*(-2*v_n*p_n*q_n/(u_n^2+v_n^2+1));
  q_n1_k= q_n+h*(-2*u_n*p_n*q_n/(u_n^2+v_n^2+1));
  
  u_n1_k1= 0;
  v_n1_k1= 0;
  p_n1_k1= 0;
  q_n1_k1= 0;
   
  while (abs(u_n1_k1-u_n1_k) >= toler && abs(v_n1_k1-v_n1_k) >= toler && abs(p_n1_k1-p_n1_k) >= toler && abs(q_n1_k1-q_n1_k) >= toler && k<10)
    #fprintf('euler atras iter de paso \n');

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

function [uOut,vOut,pOut,qOut]=paso_AB_2(u_1, v_1, p_1, q_1, u, v, p, q, h)	
	uOut= u_1 + h*((3/2)*p_1 - (1/2)*p);
	vOut= v_1 + h*((3/2)*q_1 - (1/2)*q);
	pOut= p_1 + h*((3/2)*((-2*v_1*p_1*q_1)/(u_1^2+v_1^2+1)) - (1/2)*(-2*v*p*q)/(u^2+v^2+1));
	qOut= q_1 + h*((3/2)*((-2*u_1*p_1*q_1)/(u_1^2+v_1^2+1)) - (1/2)*(-2*u*p*q)/(u^2+v^2+1));
end

function n=norma(x)
  n= (x(1)^2+x(2)^2+x(3)^2+x(4)^2)^(1/2);
end

function norma_error_max=normaMaxError(X,Y)
      norma_error_max= 0;
      for i=1:101
        error= X(i,:)-Y(i,:);
        norma_error= norma(error);
        if (norma_error>norma_error_max)
          norma_error_max= norma_error;
        end
      end 
end

function [u,v,p,q]=condIniciales()
        #u= 1; v= 3; p= 1; q= 2;
        u= 0; v= 0; p= 1; q= 2;
end

toler= 10^-4;
#r = [0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,];
#r = [0.0025, 0.005, 0.01, 0.025, 0.05]
r = [0.005, 0.01, 0.025, 0.05];

j=0;
for h=r
      fprintf('h=%d \n', h);
      
      a=0;
      b=0.1;

      RES_EA=[];
      RES_AT=[];
      RES_T=[];
      RES_H=[];
      RES_RK=[];
      
      j=j+1;
  
      k=0.05/h;
      
      fprintf('ode45 \n');
      #ode45#
      [u,v,p,q]= condIniciales();
      [odeT,odeRes] = ode45("f", (a:h:b)', [u,v,p,q]);
      #odeRes = lsode("f",[u,v,p,q],(a:h:b));
      
      fprintf('EULER HACIA ADELANTE \n');
      #EULER HACIA ADELANTE#
      [u,v,p,q]= condIniciales();
      RES_EA(1,:)=[u,v,p,q];
      
      for i=1:((b-a)/h)
        [u,v,p,q]=paso_euler_adelante(u,v,p,q,h,toler);
        RES_EA(i+1,:)=[u,v,p,q];
      end
      
      EA(j)= norma(odeRes(k,:)-RES_EA(k,:));
      fprintf('error de euler adelante en h=%d: %d \n',-log10(h),-log10(EA(j)));

      fprintf('EULER HACIA ATRAS \n');
      #EULER HACIA ATRAS#
      [u,v,p,q]= condIniciales();
      RES_AT(1,:)=[u,v,p,q];
      
      for i=1:((b-a)/h)
        [u,v,p,q]=paso_euler_atras(u,v,p,q,h,toler);
        RES_AT(i+1,:)=[u,v,p,q];
      end

      AT(j)= norma(odeRes(k,:)-RES_AT(k,:));
      fprintf('error de euler atras en h=%d: %d \n',-log10(h),-log10(AT(j)));
      
      fprintf('TRAPECIO \n');
      #TRAPECIO
      [u,v,p,q]= condIniciales();
      RES_T(1,:)=[u,v,p,q];

      for i=1:((b-a)/h)
        [u,v,p,q]=paso_trapecio(u,v,p,q,h,toler);
        RES_T(i+1,:)=[u,v,p,q];
      end
      
      T(j)= norma(odeRes(k,:)-RES_T(k,:));
      fprintf('error de trapecio en h=%d: %d \n',-log10(h),-log10(T(j)));
      
      fprintf('HEUN \n');
      #HEUN
      [u,v,p,q]= condIniciales();
      RES_H(1,:)=[u,v,p,q];

      for i=1:((b-a)/h)
        [u,v,p,q]=paso_heun(u,v,p,q,h,toler);
        RES_H(i+1,:)=[u,v,p,q];
      end
      
      H(j)= norma(odeRes(k,:)-RES_H(k,:));
      fprintf('error de heun en h=%d: %d \n',-log10(h),-log10(H(j)));
      
      fprintf('RK \n');
      #RK#
      [u,v,p,q]= condIniciales();
      RES_RK(1,:)=[u,v,p,q];

      for i=1:((b-a)/h)
        [u,v,p,q]=paso_RK(u,v,p,q,h,toler);
        RES_RK(i+1,:)=[u,v,p,q];
      end
      
      RK(j)= norma(odeRes(k,:)-RES_RK(k,:));
      fprintf('error de rk en h=%d: %d \n',-log10(h),-log10(RK(j)));
      
       
      fprintf('AB \n');
      #AB#
      [u,v,p,q]= condIniciales();
      RES_AB(1,:)=[u,v,p,q];
      [u_1,v_1,p_1,q_1]=paso_euler_adelante(u,v,p,q,h);
      RES_AB(2,:)=[u_1,v_1,p_1,q_1];
     
      for i=2:((b-a)/h)
        
        [u_2,v_2,p_2,q_2]=paso_AB_2(u_1,v_1,p_1,q_1,u,v,p,q,h);
        RES_AB(i+1,:)=[u_2,v_2,p_2,q_2];
        
        u = u_1;
        v = v_1;
        p = p_1;
        q = q_1;
        
        u_1 = u_2;
        v_1 = v_2;
        p_1 = p_2;
        q_1 = q_2;
      end
      
      AB(j)= norma(odeRes(k,:)-RES_AB(k,:));
      fprintf('error de ab en h=%d: %d \n',-log10(h),-log10(AB(j)));
end
      
#plot (r,EA,";Euler adelante;",r,AT,";Euler atras;",r,T,";Trapecio;",r,H,";Heun;",r,RK,";RK;");
#plot (r,T,";Trapecio;",r,RK,";RK;");
plot (-log10(r),-log10(EA),";Euler adelante;",-log10(r),-log10(AT),";Euler atras;",-log10(r),-log10(T),";Trapecio;",-log10(r),-log10(H),";Heun;",-log10(r),-log10(RK),";RK;",-log10(r),-log10(AB),";AB;");
