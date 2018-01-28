clc;
clear all;
ch=menu('Transmission Analysis & design','Ferranti effect','Different line models','Voltagge regulation & Efficiency','about me');
    if (ch == 1)
        disp('This illustrates the Ferranti effect');
        disp('it simulates the effect by varying the length of line');
        disp('to 5000km in steps of 10 km(by defualt)');
        disp('and plots the sending end phasor Vs');
        f  = 50;
        Vr = (220e3/sqrt(3));
        l = input('enter length of line in km   ');
        if isempty(l)
           l = 5000;
        end
        k = 1;
        r = input('enter vslue of line resistance in ohm/km   ');
        if isempty(r)
            r = 0.125;
        end
        x = input('enter the value of XL in ohm/km   ');
        if isempty(x)
            x = 0.4;
        end
        y = input('enter the value of Y(shunt) in ohm/km   ');
        if isempty(y)
            y = j*2.8e-6;
        end
        Z = r + j*x;
        gamma = sqrt(y*Z);
        alpha = real(gamma);
        beta = imag(gamma);
        k = 1;    
        for i = 0:10:l,
            Vs = (Vr/2)*exp(alpha*i)*exp(j*beta*i)+(Vr/2)*exp(-alpha*i)*exp(-j*beta*i);
            X(k) = real(Vs);
            Y(k) = imag(Vs);
            k = k + 1;
        end
        plot(X,Y)
    end
    if (ch==2)
        disp('program illustrates for different line loading')
        clear
        f = 50;
        l = 300;
        Z = 40 + j*125;
        Y = 1e-3;
        PR = 50e6/3;
        VR = 220e3/(sqrt(3));
        pfload = 0.8;
        IR = PR/(VR*pfload);
        z = Z/l;
        y = Y/l;
        disp('now we calculate sending end voltage and current');
        disp('varying line loading from 10 to 300');
        ch2=menu('choose methods','short line approximation','Nominal Pi model','Exact transmission line equation','approx of exact equations','ALL');
            if(ch2==1)
            i =1;
            for l =10:10:300
            disp('short line approximation');
            Vs_shortline(i) = VR+(z*l)*IR
            Is_shortline(i) = IR
            Spf_shortline(i)= cos(angle(Vs_shortline(i)-angle(Is_shortline(i))))
            Spower_shortline(i) = real(Vs_shortline(i)*conj(Is_shortline(i)))
            point(i) = i;
            i = i + 1;
            plot(point,abs(Vs_shortline),'r',point,abs(Is_shortline),'g',point,abs(Spf_shortline),'b',point,abs(Spower_shortline),'k')
            end
            end
            if(ch2 == 2)
                disp('Nominal pi method');
                i = 1;
                for l =10:10:300
                A = 1 + (((y*l)*(z*l))/2);
                D = A;
                B = (z*l);
                C = (y*l)*(1+((z*l)*(y*l)/4));
                Vs_nominalpi(i) = A*VR +B*IR
                Is_nominalpi(i) = C*VR +D*IR
                Spf_nominalpi(i) = cos(angle(Vs_nominalpi(i))-angle(Is_nominalpi(i)))
                Spower_nominalpi(i) = real(Vs_nominalpi(i)*conj(Is_nominalpi(i)))
                point(i) = i;
                i = i + 1;
                plot(point,abs(Vs_nominalpi),'r',point,abs(Is_nominalpi),'g',point,abs(Spf_nominalpi),'b',point,abs(Spower_nominalpi),'k')
                end
            end
                if(ch2 == 3)
                 i = 1;
                 for l = 10:10:300
                 disp('Exact trans line equation');
                 Zc = sqrt(z/y);
                 gmma = sqrt(z*y);
                 Vs_exact(i) = cosh(gmma*l)*VR +Zc*sinh(gamma*l)*IR
                 Is_exact(i) = (1/Zc)*sinh(gmma*l) + cosh(gmma*l)*IR
                 Spf_exact(i) = cos(angle(Vs_exact(i)-angle(Is_exact(i))))
                 Spower_exact(i) = real(Vs_exact(i)*conj(Is_exact(i)))
                 point(i) = i;
                 i = i + 1;
                 plot(point,abs(Vs_exact),'r',point,abs(Is_exact),'g',point,abs(Spf_exact),'b',point,abs(Spower_exact),'k')
                 end
                end
                if(ch2 == 4)
                  i = 1;
                  for l = 10:10:300
                  disp('Approx. exact equation');
                  A = 1 + ((z*l)*(y*l))/2;
                  D = A;
                  B = (z*l)*(1+((z*l)*(y*l)/6));
                  C = (y*l)*(1+((z*l)*(y*l)/6));
                  Vs_approx(i) = A*VR + B*IR;
                  Is_approx(i) = C*VR + D*IR;
                  Spf_approx(i) = cos(angle(Vs_approx(i)-angle(Is_approx(i))))
                  Spower_approx(i) = real(Vs_approx(i)*conj(Is_approx(i)))
                  point(i) = i;
                  i = i + 1;
                  plot(point,abs(Vs_approx),'r',point,abs(Is_approx),'g',point,abs(Spf_approx),'b',point,abs(Spower_approx),'k')
                  end
                end
               if(ch2 == 5)
                   i = 1;
                   for l = 10:10:300
                    disp('all line loads modelling...');
                    Vs_shortline(i) = VR+(z*l)*IR
                    Is_shortline(i) = IR
                    Spf_shortline(i)= cos(angle(Vs_shortline(i)-angle(Is_shortline(i))))
                    Spower_shortline(i) = real(Vs_shortline(i)*conj(Is_shortline(i)))
                    A = 1 + (((y*l)*(z*l))/2);
                    D = A;
                    B = (z*l);
                    C = Y*(1+((y*l)*(z*l))/4);
                    Vs_nominalpi(i) = A*VR +B*IR
                    Is_nominalpi(i) = C*VR +D*IR
                    Spf_nominalpi(i) = cos(angle(Vs_nominalpi(i)-angle(Is_nominalpi(i))))
                    Spower_nominalpi(i) = real(Vs_nominalpi(i)*conj(Is_nominalpi(i)))
                    Zc = sqrt(z/y);
                    gmma = sqrt(z*y);
                    Vs_exact(i) = cosh(gmma*l)*VR +Zc*sinh(gmma*l)*IR
                    Is_exact(i) = (1/Zc)*sinh(gmma*l) + cosh(gmma*l)*IR
                    Spf_exact(i) = cos(angle(Vs_exact(i)-angle(Is_exact(i))))
                    Spower_exact(i) = real(Vs_exact(i)*conj(Is_exact(i)))
                    A = 1 + ((z*l)*(y*l))/2;
                    D = A;
                    B = Z*(1+(z*l)*(y*l)/6);
                    C = Y*(1+(z*l)*(y*l)/6);
                    Vs_approx(i) = A*VR + B*IR;
                    Is_approx(i) = C*VR + D*IR;
                    Spf_approx(i) = cos(angle(Vs_approx(i)-angle(Is_approx(i))))
                    Spower_approx(i) = real(Vs_approx(i)*conj(Is_approx(i)))
                    point(i) = i;
                    i = i + 1;
                    plot(point,abs(Vs_shortline),'r',point,abs(Vs_nominalpi),'g',point,abs(Vs_exact),'b',point,abs(Vs_approx),'k')
                    plot(point,abs(Is_shortline),'r',point,abs(Is_nominalpi),'g',point,abs(Is_exact),'b',point,abs(Is_approx),'k')
                    plot(point,abs(Spf_shortline),'r',point,abs(Spf_nominalpi),'g',point,abs(Spf_exact),'b',point,abs(Spf_approx),'k')
                    plot(point,abs(Spower_shortline),'r',point,abs(Spower_nominalpi),'g',point,abs(Spower_exact),'b',point,abs(Spower_approx),'k')
                   end
               end
    end
    if (ch==3)
        disp('you selected Voltage regulation');
        disp('Now select pf property');
        ch3=menu('select Power factor','lagging and unity','leading')
        if(ch3==1)
            disp('you select a vol. reg. for lag or unity pf   ');
            pfr2 = input('enter the value of receiving side pf   ');
            r2 = input('enter value of resistance in ohms/km  ');
            lt = input('enter the length of trans. line   ');
            l2 = input('enter the value of inductance in mili Henry/km   ');
            x2 = 314*l2*1e-3*lt;
            r21=r2*lt;
            Vr2 = input('enter the value of supply voltage in V of receiving side   ');
            P2 = input('enter the value of power delivered to load in Watt   ');
            I2 = P2/(pfr2*Vr2);
            Vs2 = (((((I2*r21)+(Vr2*pfr2))^2)+((Vr2*sin(acos(pfr2))+(I2*x2))^2))^(1/2));
            V_reg = ((Vs2-Vr2)/Vr2) *100;
            fprintf('Voltage regulation by method1 is %0.4f\n',V_reg);
            Vreg = (((I2*r21*pfr2)+(I2*x2*sin(acos(pfr2))))/Vr2)*100;
            fprintf('Voltage regulation by another method is %0.4f',Vreg);
            disp('Tips   ...   ...   ...');
            disp('for zero voltage regulation do as below')
            disp('set tan(phi) = R/X = cot(theta)\n where phi = receiving side cosine inverse pf''');
            disp('where theta = R/Z ');
            disp('\n');
            eta = ((Vr2*I2*pfr2)/(Vr2*I2*pfr2+((I2^2)*(r21*1e-3))))*100;
            fprintf('efficiency is.....%0.2f',eta);
        end
        if(ch3==2)
            disp('you selected a vol.reg for leading pf  ');
            pfr3 = input('enter the value of pf load side '  );
            r2 = input('enter value of resistance in ohms/km  ');
            lt = input('enter the length of trans. line in km  ');
            l2 = input('enter the value of inductance in mili Henry/km   ');
            x2 = 2*pi*50*l2*1e-3*lt;
            r21=r2*lt;
            Vr2 = input('enter the value of supply voltage in V of receiving side   ');
            P2 = input('enter the value of power delivered to load in Watt   ');
            I2 = P2/(pfr2*Vr2);
            Vs2 = ((((I2*r21+Vr2*pfr3)^2)+((Vr2*sin(acos(pfr3))+(I2*x2))^2))^(1/2));
            V_reg = ((Vs2-Vr2)/Vr2) *100;
            fprintf('Voltage regulation by method1 is %0.4f\n',V_reg);
            Vreg = ((I2*r21*pfr2)-(I2*x2*sin(acos(pfr2))))/Vr2*100;
            fprintf('Voltage regulation by method2 is %0.4f',Vreg);
            disp('Tips   ...   ...   ...');
            disp('for zero voltage regulation do as below')
            disp('set tan(phi) = R/X = cot(theta)\n where phi = receiving side cosine inverse pf''');
            disp('where theta = R/Z ');
            disp('\n \n \n ');
            eta = (Vr2*I2*pfr3)/(Vr2*I2*pfr3+((I2^2)*r21))*100;
        end
		if (ch==4)
		disp('Hello i am Darsh Gajjar');
		disp('Hope you like my work...');
		disp('Thanks...................');
    end
            
            
            