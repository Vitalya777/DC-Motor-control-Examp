classdef Line2d
    properties
        p1 (1,2) {mustBeNumeric} = [0,0];
        p2 (1,2) {mustBeNumeric} = [1,1];
    end
    properties (Dependent)
        e (1,2) {mustBeNumeric}
        L (1,1) {mustBeNumeric}
        angDeg (1,1) {mustBeNumeric}
        k (1,1) {mustBeNumeric}
        b (1,1) {mustBeNumeric}
    end

    methods

        function obj = Line2d(point1,point2)
            switch nargin
                case 1
                    obj.p2 = point1;
                case 2
                    obj.p1 = point1;
                    obj.p2 = point2;
             end
        end

        function obj = evalY(obj,val)
            % находит координату X для прямой по Y
            obj = valY(obj,val);
        end

        function obj = evalX(obj,val)
            % находит координату Y для прямой по X
            obj = valX(obj,val);
        end

        function obj = LineTwoPoint(obj,point1,point2)
            % задание прямой по двум точкам
            obj.p1 = point1;
            obj.p2 = point2;
        end

        function obj = LineCoeffX(obj,c1,c2,x)
            % задание прямой по коэффициентам k, b и X
            obj.p1(1) = x;
            obj.p1(2) = c1*x+c2; 
            obj.k = c1;
            obj.b = c2;
        end

        function obj = LinePointVec(obj,point1,vec)
            % задание прямой по направляющему вектору и точке
            obj.p1 = point1;
            obj.e = vec;
        end

        function obj = LinePointLengthAngle(obj,point,l,ang,dim)
            % задание прямой по точке, длинне и углу
            n = nargin;
            obj.p1 = point;
            obj.L = l;
            if n == 5
                cnt = containers.Map({'rad','db/dec'}, {true,false});
                if cnt(dim)
                    obj.angDeg = rad2deg(ang);
                else
                    obj.angDeg = rad2deg(atan(1/ang));
                end
            else
            obj.angDeg = ang;
            end
        end

        function eqn = dispGenEqn(obj)
            % вывод уравнения для прямой
            syms x y
            eqn = vpa((x-obj.p2(1))/obj.e(1) == (y-obj.p2(2))/obj.e(2),3);
        end

        function eqn = dispTwoPoint(obj)
            % вывод уравнения для прямой по двум точкам
            syms x y
            eqn = (y-obj.p2(2))/(obj.p2(2)-obj.p1(2)) == (x-obj.p2(1))/(obj.p2(1)-obj.p1(1));
        end

        function eqn = dispCoefEqn(obj)
            % вывод уравнения для прямой в виде y = k*x+b
            syms x y
            eqn = vpa(y == obj.k*x + obj.b,3);
        end

        function point = CrossLine(obj1,obj2)
            % нахождение точки пересечения двух прямых
             point(1) = -(obj2.e(1)*obj1.e(2)*obj1.p1(1) - obj1.e(1)*obj2.e(2)*obj2.p1(1) - obj1.e(1)*obj2.e(1)*obj1.p1(2) + obj1.e(1)*obj2.e(1)*obj2.p1(2))/(obj1.e(1)*obj2.e(2) - obj2.e(1)*obj1.e(2));
             point(2) = -(obj1.e(2)*obj2.e(2)*obj1.p1(1) - obj1.e(2)*obj2.e(2)*obj2.p1(1) - obj1.e(1)*obj2.e(2)*obj1.p1(2) + obj2.e(1)*obj1.e(2)*obj2.p1(2))/(obj1.e(1)*obj2.e(2) - obj2.e(1)*obj1.e(2));
        end

        function angdb(obj)
            % нахождение угла наклона прямой в ДБ/дек
            db = evalX(obj,1) - evalX(obj,0);
            disp(string(db) + "ДБ/дек")
        end

        function plot(obj)
            % граффик прямой
            plot(obj.p1(1),obj.p1(2),'.','Color','k','MarkerSize',10)
            hold on
            plot(obj.p2(1),obj.p2(2),'.','Color','k','MarkerSize',10)
            line([obj.p1(1) obj.p2(1)],[obj.p1(2) obj.p2(2)],'LineWidth',2,'Color','b')
            grid on
        end

        function obj = plus(obj1,obj2)
            obj = Line2d;
            obj.p1 = obj1.p1 + obj2.p1;
            obj.p2 = obj1.p2 + obj2.p2;
        end

        function obj = minus(obj1,obj2)
            obj = Line2d;
            obj.p1 = obj1.p1 - obj2.p1;
            obj.p2 = obj1.p2 - obj2.p2;
        end

        function obj = set.e(obj,val)
            obj.p2 = val + obj.p1;
        end

        function val = get.e(obj)
            val = obj.p2 - obj.p1;
        end

        function obj = set.L(obj,val)
            obj.p2 = YPoint(obj,val,obj.angDeg);
        end

        function val = get.L(obj)
            val = hypot(obj.e(1),obj.e(2));
        end

        function obj = set.angDeg(obj,val)
            obj.p2 = YPoint(obj,obj.L,val);
        end
        
        function val = get.angDeg(obj)
            val = rad2deg(atan(obj.e(2)/obj.e(1)));
        end

        function obj = set.k(obj,val)
            ang = rad2deg(atan(val));
            obj.p2 = YPoint(obj,obj.L,ang);
        end

        function val = get.k(obj)
            val = tan(deg2rad(obj.angDeg));
        end

        function obj = set.b(obj,val)
            obj.p1(2) = obj.k*obj.p1(1)+val;
        end

        function val = get.b(obj)
            val = obj.k*obj.p1(1)+obj.p1(2);
        end

    end

    methods (Access=private)
        function r = YPoint(obj,l,ang)
            ang = deg2rad(ang);
            r = obj.p1+[l*(1/(tan(ang)^2 + 1))^(1/2) l*tan(ang)*(1/(tan(ang)^2 + 1))^(1/2)];
        end
        function obj = valY(obj,x)
            arguments
                obj
                x {mustBeNumeric}
            end            
            obj = (obj.e(2)*x - obj.e(2)*obj.p1(1) + obj.e(1)*obj.p1(2))/obj.e(1);
        end
        function obj = valX(obj,x)
            arguments
                obj
                x {mustBeNumeric}
            end            
            obj = (obj.e(1)*x + obj.e(2)*obj.p1(1) - obj.e(1)*obj.p1(2))/obj.e(2);
        end
    end
    
end