//                             XXXXXXX       XXXXXXX     OOOOOOOOO     LLLLLLLLLLL                    CCCCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTT   SSSSSSSSSSSSSSS 
//                             X:::::X       X:::::X   OO:::::::::OO   L:::::::::L                 CCC::::::::::::CT:::::::::::::::::::::T SS:::::::::::::::S
//                             X:::::X       X:::::X OO:::::::::::::OO L:::::::::L               CC:::::::::::::::CT:::::::::::::::::::::TS:::::SSSSSS::::::S
//                             X::::::X     X::::::XO:::::::OOO:::::::OLL:::::::LL              C:::::CCCCCCCC::::CT:::::TT:::::::TT:::::TS:::::S     SSSSSSS
//                             XXX:::::X   X:::::XXXO::::::O   O::::::O  L:::::L               C:::::C       CCCCCCTTTTTT  T:::::T  TTTTTTS:::::S            
//                                X:::::X X:::::X   O:::::O     O:::::O  L:::::L              C:::::C                      T:::::T        S:::::S            
//                                 X:::::X:::::X    O:::::O     O:::::O  L:::::L              C:::::C                      T:::::T         S::::SSSS         
//                                  X:::::::::X     O:::::O     O:::::O  L:::::L              C:::::C                      T:::::T          SS::::::SSSSS    
//                                  X:::::::::X     O:::::O     O:::::O  L:::::L              C:::::C                      T:::::T            SSS::::::::SS  
//                                 X:::::X:::::X    O:::::O     O:::::O  L:::::L              C:::::C                      T:::::T               SSSSSS::::S 
//                                X:::::X X:::::X   O:::::O     O:::::O  L:::::L              C:::::C                      T:::::T                    S:::::S
//                             XXX:::::X   X:::::XXXO::::::O   O::::::O  L:::::L         LLLLLLC:::::C       CCCCCC        T:::::T                    S:::::S
//                             X::::::X     X::::::XO:::::::OOO:::::::OLL:::::::LLLLLLLLL:::::L C:::::CCCCCCCC::::C      TT:::::::TT      SSSSSSS     S:::::S
//                             X:::::X       X:::::X OO:::::::::::::OO L::::::::::::::::::::::L  CC:::::::::::::::C      T:::::::::T      S::::::SSSSSS:::::S
//                             X:::::X       X:::::X   OO:::::::::OO   L::::::::::::::::::::::L    CCC::::::::::::C      T:::::::::T      S:::::::::::::::SS 
//                             XXXXXXX       XXXXXXX     OOOOOOOOO     LLLLLLLLLLLLLLLLLLLLLLLL       CCCCCCCCCCCCC      TTTTTTTTTTT       SSSSSSSSSSSSSSS   



////////////////////////////////////////////
//	Name        :    XOLCTS   			  //
//	Author      :    Xion			      //
//	Date        :    21/04/2021           //
//	Description :    Xion OLC Tool Set.   //
////////////////////////////////////////////



#include <olc/olcPixelGameEngine.h>
#include "Noise/FastNoiseLite.h"
#define OLC_PGEX_FONT
#include "olcPGEX/olcPGEX_Font.h"
// #define OLC_PGEX_SOUND
// #include "../olcPGEX/olcPGEX_Sound.h"
#include <sstream>
#include <random>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <unordered_set>
#if defined _WIN32
  #include <windows.h>
#else
  #include <X11/Xlib.h>
#endif

#pragma region XOLCTS



// TODO -- SubMenu



// Own definiton of dPi.
inline double dPi = 3.14159265358979323846264338327950288419716939937510582;



float Map(float value,float start1,float stop1, float start2, float stop2) 
{
    return start2 + (stop2 - start2) * ((value - start1) / (stop1 - start1));
}



float mag(olc::vf2d A, olc::vf2d B) {
        return sqrt(pow((A.x - B.x),2) + pow((A.y - B.y),2));
}



float Constrain(float value, float min_, float max_) {
    if (value > max_-1) {
        value = max_;
    }
    if (value < min_+1) {
        value = min_;
    }
    return value;
}



olc::vf2d vf2dConstrain(olc::vf2d value, olc::vf2d min_, olc::vf2d max_) {
    if (value.x > max_.x) {
        value.x = max_.x;
    }
    if (value.x < min_.x) {
        value.x = min_.x;
    }
    if (value.y > max_.y) {
        value.y = max_.y;
    }
    if (value.y < min_.y) {
        value.y = min_.y;
    }
    return value;
}


float CalcAtt(float mA,float mB,float dAB,float G) {
    return G*(mA*mB/pow(dAB,2));
}



inline float DegToRad(float x) {
    return x / 180 * dPi;
}



float RadToDeg(float x) {
    return x / dPi * 180;
}



olc::vf2d CreateVectorFromAngle(float angle) {
    return olc::vf2d(cos(DegToRad(angle)),sin(DegToRad(angle)));
}



float AngleFromPoints(olc::vf2d Point1, olc::vf2d Point2) {
    return atan2(Point1.y - Point2.y, Point1.x - Point2.x);
}



bool LineCircleIntersection(olc::vf2d A, olc::vf2d B, olc::vf2d C, float R) {
    // compute the euclidean distance between A and B
    float LAB = sqrt( pow(B.x-A.x,2)+pow(B.y-A.y,2) );

    // compute the direction vector D from A to B
    float Dx = (B.x-A.x)/LAB;
    float Dy = (B.y-A.y)/LAB;

    // the equation of the line AB is x = Dx*t + Ax, y = Dy*t + Ay with 0 <= t <= LAB.

    // compute the distance between the points A and E, where
    // E is the point of AB closest the circle center (Cx, Cy)
    float t = Dx*(C.x-A.x) + Dy*(C.y-A.y);

    // compute the coordinates of the point E
    float Ex = t*Dx+A.x;
    float Ey = t*Dy+A.y;

    // compute the euclidean distance between E and C
    float LEC = sqrt(pow(Ex-C.x,2)+pow(Ey-C.y,2));

    // test if the line intersects the circle
    if( LEC < R ) {
        // compute distance from t to circle intersection point
        float dt = sqrt( pow(R,2) - pow(LEC,2));

        // compute first intersection point
        float Fx = (t-dt)*Dx + A.x;
        float Fy = (t-dt)*Dy + A.y;

        // compute second intersection point
        float Gx = (t+dt)*Dx + A.x;
        float Gy = (t+dt)*Dy + A.y;
        return true;
    } else if( LEC == R ) { // else test if the line is tangent to circle
        return false;
    }else {
        return false;
    }
    return false;
}



float RandomFloat(float MinVal, float MaxVal) {
        return MinVal + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(MaxVal-MinVal)));
    };



olc::vf2d RandomPointInCircle(olc::vf2d vLocation, float fR, float fMinR ) {
    float r = fR * sqrt(RandomFloat(Constrain(fMinR/fR,0,1),1));
    float theta = RandomFloat(0,1) * 2 * dPi;
    return olc::vf2d(vLocation.x + r * cos(theta), vLocation.y + r * sin(theta));
}



olc::vf2d RandomPointInCircumference(olc::vf2d vLocation, float fR) {
    float theta = RandomFloat(0,1) * 2 * dPi;
    return olc::vf2d(vLocation.x + fR * cos(theta), vLocation.y + fR * sin(theta));
}



float Lerp(float A, float B, float C) {
    return (1 - C) * A + C * B;
}



bool CircleCircleColision(olc::vf2d A, float fAR,olc::vf2d B, float fBR) {

    float dx = A.x - B.x;
    float dy = A.y - B.y;
    float distance = sqrt(dx * dx + dy * dy);

    if (distance < fAR + fBR) {
        return true;
    }
    return false;
}



float GetBigest(std::vector<float> Values) {
    float record = 0;
    float recordIndex = 0;
    for(int i = 0; i < Values.size(); i++) {
        if(Values[i] > record) {
            record = Values[i]; 
            recordIndex = i;
        } 
    }
    return recordIndex;
}



olc::vf2d EndFromAngle(olc::vf2d A, float Angle, float length) {
    olc::vf2d AngleV = CreateVectorFromAngle(Angle);
    return olc::vf2d((float)A.x+AngleV.x*length,(float)A.y+AngleV.y*length);
}



olc::vf2d EndFromAngle(olc::vf2d A, olc::vf2d Angle, float length) {
        return olc::vf2d((float)A.x+Angle.x*length,(float)A.y+Angle.y*length);
    }



bool LineInterserct(olc::vf2d StartA,olc::vf2d EndA,olc::vf2d StartB,olc::vf2d EndB, olc::vf2d *Output) {
        float x1 = StartA.x;
        float y1 = StartA.y;
        float x2 = EndA.x;
        float y2 = EndA.y;
        float x3 = StartB.x;
        float y3 = StartB.y;
        float x4 = EndB.x;
        float y4 = EndB.y;
        float den = (x1 -x2) * (y3 - y4) - (y1 -y2) * (x3 - x4);
        if (den == 0) {
            return false;
        }
        float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
        float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / den;
        if (t > 0 and t < 1 and u > 0) {
            float ptx = x1 + t * (x2 - x1);
            float pty = y1 + t * (y2 - y1);
            olc::vf2d pt = olc::vf2d(ptx,pty);
            if (mag(olc::vf2d(ptx,pty),olc::vf2d(x2,y2)) < mag(olc::vf2d(x1,y1),olc::vf2d(x2,y2))) {
                if (mag(olc::vf2d(ptx,pty),olc::vf2d(x3,y3)) <= mag(olc::vf2d(x3,y3),olc::vf2d(x4,y4))) {
                    Output->x = pt.x;
                    Output->y = pt.y;
                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        } else {
            return false;
        }
    }



bool InBounds(olc::vf2d TopLeft, olc::vf2d BottomRight, olc::vf2d Point) {
        if(Point.x >= TopLeft.x && Point.x <= BottomRight.x && Point.y >= TopLeft.y && Point.y <= BottomRight.y) {
            return true;
        } 
        return false;
    }



namespace XOLCTS {



    #pragma region Variables



    class Renderable;



    struct KeyboardShortcut {
        std::vector<olc::Key> vKeys;
        void (*vdActionFunction)();
        bool bHoldable = false;

        void CheckKeypresses(olc::PixelGameEngine *pge, int iKeyPressedCount) {
            if(iKeyPressedCount != vKeys.size()) return;
            for(auto &i : vKeys) {
                if(!pge->GetKey(i).bHeld) return;
            }
            if(!pge->GetKey(vKeys[vKeys.size()-1]).bPressed && !bHoldable) return;
            vdActionFunction();
        }
    };



    std::vector<KeyboardShortcut*> vKeyboardShortcuts;



    std::vector<Renderable*> vcRenderableList;



    struct MoveRequest {
        int iTime;
        olc::vf2d vOrigin;
        olc::vf2d vDestination;
        olc::vf2d vIncrement;
        int iID;
        XOLCTS::Renderable* e;
    };



    std::vector<MoveRequest*> vcMoveRequestList;



    struct TimeMesure {
        std::string Name;
		std::chrono::steady_clock::time_point tpStartPoint = std::chrono::steady_clock::now();
		std::chrono::steady_clock::time_point tpEndPoint = std::chrono::steady_clock::now();
        double dTime = 0;
        int iDivisor = 0;
    };



    std::vector<TimeMesure*> vTimeMesurements;



    // Spline Struct
    struct Spline {
        std::array<olc::vf2d, 4> vPoints;

        olc::vf2d GetSplinePoint(float t) {
            int p0, p1, p2,p3;
            p1 = (int)t + 1;
            p2 = p1 + 1;
            p3 = p2 + 1;
            p0 = p1 - 1;

            float tt = t * t;
            float ttt = tt * t;

            float q1 = -ttt + 2 * tt - t;
            float q2 = 3*ttt - 5*tt + 2;
            float q3 = -3*ttt + 4*tt + t;
            float q4 = ttt - tt;

            float tx =  vPoints[p0].x * q1 + 
                        vPoints[p1].x * q2 +
                        vPoints[p2].x * q3 + 
                        vPoints[p3].x * q4;

            float ty =  vPoints[p0].y * q1 + 
                        vPoints[p1].y * q2 +
                        vPoints[p2].y * q3 + 
                        vPoints[p3].y * q4;

            return olc::vf2d(tx*0.5,ty*0.5);
        }
    
        olc::vf2d DrawSpline(float Precision, olc::PixelGameEngine* pge) {
            olc::vf2d vLastPos = vPoints[1];
            std::vector<olc::vf2d> vPointsRecord;
            for(float t = 0; t < 1; t+= Precision) {
                olc::vf2d v = GetSplinePoint(t);
                vPointsRecord.push_back(v);
                pge->DrawLine(vLastPos,v);
                vLastPos = v;

            }
            pge->DrawLine(vPointsRecord[vPointsRecord.size()-1],vPoints[2]);
            return vPointsRecord[(int)((vPointsRecord.size()-1)/3*2)];
        }

        olc::vf2d DrawSpline(float Precision, olc::PixelGameEngine* pge, olc::Pixel p) {
            olc::vf2d vLastPos = vPoints[1];
            std::vector<olc::vf2d> vPointsRecord;
            for(float t = 0; t < 1; t+= Precision) {
                olc::vf2d v = GetSplinePoint(t);
                vPointsRecord.push_back(v);
                pge->DrawLine(vLastPos,v, p);
                vLastPos = v;

            }
            pge->DrawLine(vPointsRecord[vPointsRecord.size()-1],vPoints[2],p);
            return vPointsRecord[(int)((vPointsRecord.size()-1)/3*2)];
        }

        olc::vf2d DrawSpline(float Precision, olc::PixelGameEngine* pge, olc::Pixel p, olc::Pixel p2);

        olc::vf2d DrawArrowedSpline(float Precision, olc::PixelGameEngine* pge, olc::Pixel p, int offset) {
            olc::vf2d vLastPos = vPoints[1];
            olc::vf2d vActualPos = {0,0};
            std::vector<olc::vf2d> vPointsRecord;
            for(float t = 0; t < 1; t+= Precision) {
                vActualPos = GetSplinePoint(t);
                vPointsRecord.push_back(vActualPos);
                pge->DrawLine(vLastPos,vActualPos, p);
                vLastPos = vActualPos;
            }
            
            olc::vf2d vMiddleStart = vPointsRecord[(int)((vPointsRecord.size()-1)/3*2)];
            olc::vf2d vMiddleEnd = vPointsRecord[(int)((vPointsRecord.size()-1)/3*2+1)];
            olc::vf2d OffsetEnd = EndFromAngle(vMiddleStart,180+RadToDeg(AngleFromPoints(vMiddleStart,vMiddleEnd)),mag(vMiddleEnd,vMiddleStart));
            pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vMiddleStart,vMiddleEnd))+30,5),p);
            pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vMiddleStart,vMiddleEnd))-30,5),p);
            pge->DrawLine(vPointsRecord[vPointsRecord.size()-1],vPoints[2],p);
            return vPointsRecord[(int)((vPointsRecord.size()-1)/3*2)];
        }

        olc::vf2d DrawArrowedSpline(float Precision, olc::PixelGameEngine* pge, olc::Pixel p, olc::Pixel p2, int offset);


    };



    // XOLCTS Menu Class Vars
    int iMenuOpened = -1;



    // XOLCTS SubMenu ClassVar
    int iSubMenuOpened = -1;
    


    // XOLCTS Menu Class Vars
    bool bTopMenuShown = true;



    // Prints The name of the time mesure when it starts and ends
    bool bTimeMesureDebug = false;



    // Basicly a draw Offset
    olc::vf2d vShake = {0,0};



    // Basicly a draw Scale offset
    olc::vf2d vZoom = {1,1};


    // MaxThreadCount
    static unsigned long long thread_nr = 0;



    // MaxNeededThreadCount
    static unsigned long long threadSearch_nr = 100;

    
    
    // Thread for GetMaxThreadCount
    pthread_mutex_t mutex_;



    // Variable to check if an element is actualy in use
    bool bUsesRenderable = false;



    // Default XOLCTS color for background
    olc::Pixel pBackgroundColor = olc::Pixel(38,37,37,255);



    // Used for Mouse Delta
    olc::vf2d vLastMousePos = olc::vf2d(0,0);



    // Thing for Linux clipboard 
    int XA_STRING = 31;



    #pragma endregion Variables




    #pragma region Functions



    // Clear Console
    void Clear() {
        #if defined _WIN32
            system("cls");
            //clrscr(); // including header file : conio.h
        #elif defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
            // system("clear");
            std::cout<< u8"\033[2J\033[1;1H"; //Using ANSI Escape Sequences 
        #elif defined (__APPLE__)
            system("clear");
        #endif
        return;
    }



    // Set Cursor Position
    void setCursorPosition(int x, int y) {
        #if defined _WIN32
            static const HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
            std::cout.flush();
            COORD coord = { (SHORT)x, (SHORT)y };
            SetConsoleCursorPosition(hOut, coord);
        #elif defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
            printf("\033[%d;%dH",x+1,y+1);
        #endif
    }



    // Draws a simple Arrow
    void DrawArrow(olc::vf2d vStart, olc::vf2d vEnd, olc::PixelGameEngine* pge, int offset) {
        // olc::vf2d OffsetEnd = vEnd-(vEnd.norm()*(offset+4));
        olc::vf2d OffsetEnd = EndFromAngle(vStart,180+RadToDeg(AngleFromPoints(vStart,vEnd)),mag(vStart,vEnd)-offset);
        pge->DrawLine(vStart,OffsetEnd);
        pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vStart,vEnd))+30,5));
        pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vStart,vEnd))-30,5));
        return;
    }



    // Draws a simple Arrow
    void DrawArrow(olc::vf2d vStart, olc::vf2d vEnd, olc::PixelGameEngine* pge, int offset, olc::Pixel p) {
        // olc::vf2d OffsetEnd = vEnd-(vEnd.norm()*(offset+4));
        olc::vf2d OffsetEnd = EndFromAngle(vStart,180+RadToDeg(AngleFromPoints(vStart,vEnd)),mag(vStart,vEnd)-offset);
        pge->DrawLine(vStart,OffsetEnd,p);
        pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vStart,vEnd))+30,5),p);
        pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vStart,vEnd))-30,5),p);
        return;
    }



    // Clears the time mesurements Vector
    void ClearTimeMesures() {
        vTimeMesurements.clear();
        return;
    }



    // Starts a new time mesure
    void StartTimeMesure(std::string sMesureName) {
        if(bTimeMesureDebug) std::cout << "Starting : " + sMesureName<< std::endl;
        for(auto &i : vTimeMesurements) {
            if(i->Name == sMesureName) {
                i->tpStartPoint = std::chrono::steady_clock::now();
                return;
            }
        }
		vTimeMesurements.push_back(new TimeMesure{sMesureName, std::chrono::steady_clock::now()});
        return;
    }

    

    // Stops time mesure and adds the result to the time mesure vector 
    void StopTimeMesure(std::string sMesureName) {
        if(bTimeMesureDebug) std::cout << "Ending : " + sMesureName<< std::endl;
        for(auto &i : vTimeMesurements) {
            if(i->Name == sMesureName) {
                i->iDivisor++;
		        i->tpEndPoint = std::chrono::steady_clock::now();
                i->dTime += std::chrono::duration_cast<std::chrono::microseconds>(i->tpEndPoint - i->tpStartPoint).count();
                return;
            }
        }
        return;
    }



    // Prints the time mesures
    void PrintTimeMesures() {
        setCursorPosition(0,0);
        int iMaxSize = 0;
        Clear();
        for(auto &i : vTimeMesurements) {
            if(i->Name.length() > iMaxSize) {
                iMaxSize = i->Name.length();
            }
        }

        // std::cout << iMaxSize << std::endl;
        std::vector<std::string> vNamesOfTimeMesures;

        for(auto &i : vTimeMesurements) {
            std::string sTreatedName = i->Name;
            while(sTreatedName.length() < iMaxSize+4) {
                sTreatedName = sTreatedName+" ";
            }
            sTreatedName += ":   ";
            vNamesOfTimeMesures.push_back(sTreatedName);
        }


        std::cout << "----------------Time Mesures----------------" << std::endl;
        for(int i = 0; i < vNamesOfTimeMesures.size(); i++) {
            TimeMesure* tmCurrentTimeMesure = vTimeMesurements[i];
            std::cout << vNamesOfTimeMesures[i]+std::to_string(tmCurrentTimeMesure->dTime/tmCurrentTimeMesure->iDivisor)+"µs ("+std::to_string(tmCurrentTimeMesure->iDivisor)+" as divisor)" << std::endl;
        }
        std::cout << "--------------------------------------------" << std::endl;
        return;
    }



    //This function checks if the Point is in bound of the defined rectangle.
    bool InBounds(olc::vf2d TopLeft, olc::vf2d BottomRight, olc::vf2d Point) {

        return (Point.x >= TopLeft.x && Point.x <= BottomRight.x && Point.y >= TopLeft.y && Point.y <= BottomRight.y);
    
    };



    // // This function is used to map values from a certain range to another.
    // float Map(float value,float start1,float stop1, float start2, float stop2) {
    //     return start2 + (stop2 - start2) * ((value - start1) / (stop1 - start1));
    // }



    // Advanced function for getting noise using the FastNoise lib.
    float GetNoiseAdvanced(FastNoiseLite NoiseObject, float x,float y,float z, float fFrequency, float octaves, float persistence, float lacunarity) {
        float Value = 0;
        float Amplitude = 1;
        float TotalAmplitude = 0;
        float Frequency = fFrequency;
        for(int i = 0; i < octaves; i++) {
            Value += NoiseObject.GetNoise(x*Frequency,y*Frequency,z*Frequency)*Amplitude;
            TotalAmplitude += Amplitude;
            Amplitude *= persistence;
            Frequency *= lacunarity;
        }
        return Value/TotalAmplitude;
    }



    // Lerps between two values
    float Lerp(float A, float B, float C) {
	    return (1 - C) * A + C * B;
    }



    // Lerps between two colors.
    olc::Pixel ColorLerp(olc::Pixel a, olc::Pixel b, float C) {
        return(olc::Pixel(Lerp(a.r,b.r,C),Lerp(a.g,b.g,C),Lerp(a.b,b.b,C)));
    }



    // Generates a random Char
    std::string RandomString(int max_length) {
        std::string possible_characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        std::random_device rd;
        std::mt19937 engine(rd());
        std::uniform_int_distribution<> dist(0, possible_characters.size()-1);
        std::string ret = "";
        for(int i = 0; i < max_length; i++){
            int random_index = dist(engine); //get index between 0 and possible_characters.size()-1
            ret += possible_characters[random_index];
        }
        return ret;
    }



    // Painfull key to key name
    std::string GetKeyName(olc::Key k) {
        switch(k) {
            default:
                return "Unknown";
                break;
            case olc::Key::A :
                return "A";
                break;
            case olc::Key::B :
                return "B";
                break;
            case olc::Key::BACK :
                return "Back";
                break;
            case olc::Key::C :
                return "C";
                break;
            case olc::Key::CAPS_LOCK :
                return "Caps Lock";
                break;
            case olc::Key::COMMA :
                return "Comma";
                break;
            case olc::Key::CTRL :
                return "Ctrl";
                break;
            case olc::Key::D :
                return "D";
                break;
            case olc::Key::DOWN :
                return "Down";
                break;
            case olc::Key::E :
                return "E";
                break;
            case olc::Key::END :
                return "End";
                break;
            case olc::Key::ENTER :
                return "Enter";
                break;
            case olc::Key::EQUALS :
                return "=";
                break;
            case olc::Key::ESCAPE :
                return "Escape";
                break;
            case olc::Key::F10 :
                return "F10";
                break;
            case olc::Key::F11 :
                return "F11";
                break;
            case olc::Key::F12 :
                return "F12";
                break;
            case olc::Key::F1 :
                return "F1";
                break;
            case olc::Key::F2 :
                return "F2";
                break;
            case olc::Key::F3 :
                return "F3";
                break;
            case olc::Key::F4 :
                return "F4";
                break;
            case olc::Key::F5 :
                return "F5";
                break;
            case olc::Key::F6 :
                return "F6";
                break;
            case olc::Key::F7 :
                return "F7";
                break;
            case olc::Key::F8 :
                return "F8";
                break;
            case olc::Key::F9 :
                return "F9";
                break;
            case olc::Key::F :
                return "F";
                break;
            case olc::Key::G :
                return "G";
                break;
            case olc::Key::H :
                return "H";
                break;
            case olc::Key::HOME :
                return "Home";
                break;
            case olc::Key::I :
                return "I";
                break;
            case olc::Key::INS :
                return "Insert";
                break;
            case olc::Key::J :
                return "J";
                break;
            case olc::Key::K:
                return "K";
                break;
            case olc::Key::L :
                return "L";
                break;
            case olc::Key::LEFT :
                return "Left";
                break;
            case olc::Key::M :
                return "M";
                break;
            case olc::Key::MINUS :
                return "-";
                break;
            case olc::Key::N :
                return "N";
                break;
            case olc::Key::NP0 :
                return "NumPad0";
                break;
            case olc::Key::NP1 :
                return "NumPad1";
                break;
            case olc::Key::NP2 :
                return "NumPad2";
                break;
            case olc::Key::NP3 :
                return "NumPad3";
                break;
            case olc::Key::NP4 :
                return "NumPad4";
                break;
            case olc::Key::NP5 :
                return "NumPad5";
                break;
            case olc::Key::NP6 :
                return "NumPad6";
                break;
            case olc::Key::NP7 :
                return "NumPad7";
                break;
            case olc::Key::NP8 :
                return "NumPad8";
                break;
            case olc::Key::NP9 :
                return "NumPad9";
                break;
            case olc::Key::NP_ADD :
                return "NumPad_Add";
                break;
            case olc::Key::NP_DECIMAL :
                return "NumPad_Decimal";
                break;
            case olc::Key::NP_DIV :
                return "NumPad_Div";
                break;
            case olc::Key::NP_MUL :
                return "NumPad_Mul";
                break;
            case olc::Key::NP_SUB :
                return "NumPad_Sub";
                break;
            case olc::Key::O :
                return "O";
                break;
            case olc::Key::P :
                return "P";
                break;
            case olc::Key::PAUSE :
                return "Pause";
                break;
            case olc::Key::PERIOD :
                return ".";
                break;
            case olc::Key::PGDN :
                return "Page_Down";
                break;
            case olc::Key::PGUP :
                return "Page_Up";
                break;
            case olc::Key::Q :
                return "Q";
                break;
            case olc::Key::R :
                return "R";
                break;
            case olc::Key::RETURN:
                return "Return";
                break;
            case olc::Key::RIGHT :
                return "Right";
                break;
            case olc::Key::S :
                return "S";
                break;
            case olc::Key::SCROLL:
                return "Scroll";
                break;
            case olc::Key::SHIFT :
                return "Shift";
                break;
            case olc::Key::SPACE :
                return "Space";
                break;
            case olc::Key::T :
                return "T";
                break;
            case olc::Key::TAB :
                return "Tab";
                break;
            case olc::Key::U :
                return "U";
                break;
            case olc::Key::UP :
                return "Up";
                break;
            case olc::Key::V :
                return "V";
                break;
            case olc::Key::W :
                return "W";
                break;
            case olc::Key::X :
                return "X";
                break;
            case olc::Key::Y :
                return "Y";
                break;
            case olc::Key::Z :
                return "Z";
                break;
        }
        return "Error";
    }



    // Checks for keyboard shortcuts and execute their action if needed.
    void CheckKeyboardShortcuts(olc::PixelGameEngine *pge) {
        if(bUsesRenderable) return;
        int iKeypressCount = 0;
        bool bAdded[olc::mapKeys.size()];
        for(std::map<size_t, uint8_t>::iterator it = olc::mapKeys.begin(); it != olc::mapKeys.end(); ++it) {
            if(pge->GetKey(olc::Key(it->second)).bHeld && !bAdded[olc::Key(it->second)]) {bAdded[olc::Key(it->second)] = true; iKeypressCount++;}
        }

        for(auto &i : vKeyboardShortcuts) {
            i->CheckKeypresses(pge,iKeypressCount);
        }
    }



    // Adds a keyboard shortcut
    KeyboardShortcut* AddKeyboardShortcut(std::vector<olc::Key> vKeys, void (*vdActionFunction)()) {
        KeyboardShortcut* i = new KeyboardShortcut{vKeys,vdActionFunction};
        vKeyboardShortcuts.push_back(i);
        return i;
    }



    // Gets the angle between <Start> and <End> as a vector.
    olc::vf2d AngleVectorFromPoints(olc::vf2d Start, olc::vf2d End) {
        return CreateVectorFromAngle(180+RadToDeg(AngleFromPoints(Start,End)));
    }



    // Just rounds a vector bro... do i even explain what rounding is?
    olc::vf2d RoundVector(olc::vf2d v) {
        return olc::vf2d(round(v.x),round(v.y));
    }



    // Useless Function. Used by the GetMaxThreadCount Function
    void* inc_thread_nr(void* arg) {
        (void*)arg;
        pthread_mutex_lock(&mutex_);
        thread_nr ++;
        pthread_mutex_unlock(&mutex_);
        // sleep(300000);
		return nullptr;
    } 



    // Returns the max number of threads the processor can handle
    int GetMaxThreadCount() {
		int err;
		int cnt = 0;
		pthread_t pid[threadSearch_nr+1];
		pthread_mutex_init(&mutex_, NULL);
		while (cnt < threadSearch_nr) {
			err = pthread_create(&pid[cnt], NULL, inc_thread_nr, NULL);
			if (err != 0) {
				break;
			}
		    pthread_join(pid[cnt], NULL);
			cnt++;
		}
		pthread_mutex_destroy(&mutex_);

        return thread_nr;
    }



    // Splits A String
    std::vector<std::string> SplitString(std::string s, std::string delimiter) {
        size_t pos_start = 0, pos_end, delim_len = delimiter.length();
        std::string token;
        std::vector<std::string> res;

        while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
            token = s.substr (pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            res.push_back (token);
        }
        res.push_back (s.substr (pos_start));
        return res;
    }



    std::string ReplaceString(std::string s, char OriginalChar, char ReplaceChar) {
        for(int i = 0; i < s.size(); i++) {
            if(s[i] == OriginalChar) s[i] = ReplaceChar;
        }
        return s;
    }



    // Draws a Spline
    olc::vf2d DrawSpline(std::array<olc::vf2d,4> Points, float Precision, olc::PixelGameEngine* pge) {
        Spline s;
        s.vPoints = Points;
        return s.DrawSpline(Precision,pge);
    }



    // Draws a Spline
    olc::vf2d DrawSpline(std::array<olc::vf2d,4> Points, float Precision, olc::PixelGameEngine* pge, olc::Pixel p) {
        Spline s;
        s.vPoints = Points;
        return s.DrawSpline(Precision,pge,p);
    }
    
    
    
    // Draws a Spline
    olc::vf2d DrawSpline(std::array<olc::vf2d,4> Points, float Precision, olc::PixelGameEngine* pge, olc::Pixel p, olc::Pixel p2) {
        Spline s;
        s.vPoints = Points;
        return s.DrawSpline(Precision,pge,p,p2);
    }
    
    

    // Draws a Spline
    olc::vf2d DrawArrowedSpline(std::array<olc::vf2d,4> Points, float Precision, olc::PixelGameEngine* pge, olc::Pixel p, int offset) {
        Spline s;
        s.vPoints = Points;
        return s.DrawArrowedSpline(Precision,pge,p, offset);
    }



    // Draws a Spline
    olc::vf2d DrawArrowedSpline(std::array<olc::vf2d,4> Points, float Precision, olc::PixelGameEngine* pge, olc::Pixel p,olc::Pixel p2,int offset) {
        Spline s;
        s.vPoints = Points;
        return s.DrawArrowedSpline(Precision,pge,p,p2, offset);
    }


    
    // Checks if string is numeric
    bool IsNumber(const std::string& s) {
        std::string::const_iterator it = s.begin();
        while (it != s.end() && std::isdigit(*it)) ++it;
        return !s.empty() && it == s.end();
    }



    // Color from string 
    olc::Pixel ColorFromString(std::string st) {
        olc::Pixel p;
        std::vector<std::string> s = SplitString(ReplaceString(ReplaceString(st,'(',' '),')',' '),",");
        if(s.size() == 3) {
            if(IsNumber(s[0]) && IsNumber(s[1]) && IsNumber(s[2])) {
                p = olc::Pixel(stoi(s[0]),stoi(s[1]),stoi(s[2]));
            } else p = olc::WHITE;
        } else p = olc::WHITE;
        return p;
    }



    // HSL color to RGB
    olc::Pixel HSLtoRGB(float H, float S, float L) {
        H = H - (360*(int)(H/360));
        float C = (1 - std::fabs(2 * L - 1)) * S; // Chroma
        float HPrime = H / 60; // H'
        float X = C * (1 - std::fabs(std::fmod(HPrime, 2.f) - 1));
        float M = L - C / 2;

        float R = 0.f;
        float G = 0.f;
        float B = 0.f;

        switch (static_cast<int>(HPrime))
        {
        default:	break;
        case 0: R = C; G = X;        break; // [0, 1)
        case 1: R = X; G = C;        break; // [1, 2)
        case 2:        G = C; B = X; break; // [2, 3)
        case 3:        G = X; B = C; break; // [3, 4)
        case 4: R = X;        B = C; break; // [4, 5)
        case 5: R = C;        B = X; break; // [5, 6)
        }

        R += M;
        G += M;
        B += M;

        olc::Pixel color;
        color.r = static_cast<int>(std::round(R * 255));
        color.g = static_cast<int>(std::round(G * 255));
        color.b = static_cast<int>(std::round(B * 255));

        return color;
    }



    // RGB color to HSL
    std::array<float,3> RGBToHSL(olc::Pixel p) {
        float H,S,L;

        float r = (p.r / 255.0f);
        float g = (p.g / 255.0f);
        float b = (p.b / 255.0f);

        float min = std::min(std::min(r, g), b);
        float max = std::max(std::max(r, g), b);
        float delta = max - min;

        L = (max + min) / 2;

        if (delta == 0)
        {
            H = 0;
            S = 0.0f;
        }
        else
        {
            S = (L <= 0.5) ? (delta / (max + min)) : (delta / (2 - max - min));

            float hue;

            if (r == max)
            {
                hue = ((g - b) / 6) / delta;
            }
            else if (g == max)
            {
                hue = (1.0f / 3) + ((b - r) / 6) / delta;
            }
            else
            {
                hue = (2.0f / 3) + ((r - g) / 6) / delta;
            }

            if (hue < 0)
                hue += 1;
            if (hue > 1)
                hue -= 1;

            H = (int)(hue * 360);
        }

        return std::array<float,3>{H,S,L};
    }



    // Draw Spline Definition
    olc::vf2d Spline::DrawSpline(float Precision, olc::PixelGameEngine* pge, olc::Pixel p, olc::Pixel p2) {
        olc::vf2d vLastPos = vPoints[1];
        std::vector<olc::vf2d> vPointsRecord;
        for(float t = 0; t < 1; t+= Precision) {
            olc::vf2d v = GetSplinePoint(t);
            vPointsRecord.push_back(v);
            pge->DrawLine(vLastPos,v, ColorLerp(p,p2,t));
            vLastPos = v;

        }
        return vPointsRecord[(int)((vPointsRecord.size()-1)/3*2)];
    }


    // Draw Arrowed Spline Definition
    olc::vf2d Spline::DrawArrowedSpline(float Precision, olc::PixelGameEngine* pge, olc::Pixel p, olc::Pixel p2, int offset) {
        olc::vf2d vLastPos = vPoints[1];
        olc::vf2d vActualPos = {0,0};
        std::vector<olc::vf2d> vPointsRecord;
        for(float t = 0; t < 1; t+= Precision) {
            vActualPos = GetSplinePoint(t);
            vPointsRecord.push_back(vActualPos);
            pge->DrawLine(vLastPos,vActualPos, ColorLerp(p,p2,t));
            vLastPos = vActualPos;
        }
        
        olc::vf2d vMiddleStart = vPointsRecord[(int)((vPointsRecord.size()-1)/3*2)];
        olc::vf2d vMiddleEnd = vPointsRecord[(int)((vPointsRecord.size()-1)/3*2+1)];
        olc::vf2d OffsetEnd = EndFromAngle(vMiddleStart,180+RadToDeg(AngleFromPoints(vMiddleStart,vMiddleEnd)),mag(vMiddleEnd,vMiddleStart));
        pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vMiddleStart,vMiddleEnd))+30,5),ColorLerp(p,p2,0.66));
        pge->DrawLine(OffsetEnd,EndFromAngle(OffsetEnd,RadToDeg(AngleFromPoints(vMiddleStart,vMiddleEnd))-30,5),ColorLerp(p,p2,0.66));
        pge->DrawLine(vPointsRecord[vPointsRecord.size()-1],vPoints[2],p2);
        return vPointsRecord[(int)((vPointsRecord.size()-1)/3*2)];
    }



    // Updates the mouse delta function
    void UpdateMouseDelta(olc::vf2d vMouseLoc) {
        vLastMousePos = vMouseLoc;
    }



    //  Returns last frame Mouse pos - Current Pos
    olc::vf2d GetMouseDelta(olc::vf2d vMouseLoc) {
        return vMouseLoc - vLastMousePos;
    }



    // Gets Screen resolution
    void GetScreenResolution(int &width, int &height) {
        #if defined _WIN32
            width  = (int) GetSystemMetrics(SM_CXSCREEN);
            height = (int) GetSystemMetrics(SM_CYSCREEN);
        #else
            X11::Display* disp = X11::XOpenDisplay(NULL);
            X11::Screen*  scrn = X11::XDefaultScreenOfDisplay(disp);
            width  = scrn->width;
            height = scrn->height;
            X11::XCloseDisplay(disp);
        #endif
        }



    // Parse String To Vector
    olc::vf2d vf2dFromString(std::string s) {
        std::vector<std::string> t = SplitString(ReplaceString(ReplaceString(s,'(',' '),')',' '),",");
        return olc::vf2d(std::stof(t[0]),std::stof(t[1]));
    }



    // Check If dir Exists
    int DirExists(const char* const path)
    {
        struct stat info;

        int statRC = stat( path, &info );
        if( statRC != 0 )
        {
            if (errno == ENOENT)  { return 0; } // something along the path does not exist
            if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
            return -1;
        }

        return ( info.st_mode & S_IFDIR) ? 1 : 0;
    }



    // Converts char to lowercase
    char ASCIIToLower(char in) {
        if (in <= 'Z' && in >= 'A')
            return in - ('Z' - 'z');
        return in;
    }



    // Get text of key
    std::string GetKey(olc::Key k) {
        std::string s = "";
        switch(k) {
            default:
                break;
            case olc::Key::A :
                s = "A";
                break;
            case olc::Key::B :
                s = "B";
                break;
            case olc::Key::C :
                s = "C";
                break;
            case olc::Key::COMMA :
                s = ",";
                break;
            case olc::Key::D :
                s = "D";
                break;
            case olc::Key::E :
                s = "E";
                break;
            case olc::Key::EQUALS :
                s = "=";
                break;
            case olc::Key::F :
                s = "F";
                break;
            case olc::Key::G :
                s = "G";
                break;
            case olc::Key::H :
                s = "H";
                break;
            case olc::Key::I :
                s = "I";
                break;
            case olc::Key::J :
                s = "J";
                break;
            case olc::Key::K:
                s = "K";
                break;
            case olc::Key::L :
                s = "L";
                break;
            case olc::Key::M :
                s = "M";
                break;
            case olc::Key::MINUS :
                s = "-";
                break;
            case olc::Key::N :
                s = "N";
                break;
            case olc::Key::NP0 :
                s = "0";
                break;
            case olc::Key::NP1 :
                s = "1";
                break;
            case olc::Key::NP2 :
                s = "2";
                break;
            case olc::Key::NP3 :
                s = "3";
                break;
            case olc::Key::NP4 :
                s = "4";
                break;
            case olc::Key::NP5 :
                s = "5";
                break;
            case olc::Key::NP6 :
                s = "6";
                break;
            case olc::Key::NP7 :
                s = "7";
                break;
            case olc::Key::NP8 :
                s = "8";
                break;
            case olc::Key::NP9 :
                s = "9";
                break;
            case olc::Key::NP_ADD :
                s = "+";
                break;
            case olc::Key::NP_DECIMAL :
                s = ".";
                break;
            case olc::Key::NP_DIV :
                s = "/";
                break;
            case olc::Key::NP_MUL :
                s = "*";
                break;
            case olc::Key::NP_SUB :
                s = "-";
                break;
            case olc::Key::O :
                s = "O";
                break;
            case olc::Key::P :
                s = "P";
                break;
            case olc::Key::PERIOD :
                s = ".";
                break;
            case olc::Key::Q :
                s = "Q";
                break;
            case olc::Key::R :
                s = "R";
                break;
            case olc::Key::S :
                s = "S";
                break;
            case olc::Key::SPACE :
                s = " ";
                break;
            case olc::Key::T :
                s = "T";
                break;
            case olc::Key::U :
                s = "U";
                break;
            case olc::Key::V :
                s = "V";
                break;
            case olc::Key::W :
                s = "W";
                break;
            case olc::Key::X :
                s = "X";
                break;
            case olc::Key::Y :
                s = "Y";
                break;
            case olc::Key::Z :
                s = "Z";
                break;
            case olc::Key(27) :
                s = "0";
                break;
            case olc::Key(28) :
                s = "1";
                break;
            case olc::Key(29) :
                s = "2";
                break;
            case olc::Key(30) :
                s = "3";
                break;
            case olc::Key(31) :
                s = "4";
                break;
            case olc::Key(32) :
                s = "5";
                break;
            case olc::Key(33) :
                s = "6";
                break;
            case olc::Key(34) :
                s = "7";
                break;
            case olc::Key(35) :
                s = "8";
                break;
            case olc::Key(36) :
                s = "9";
                break;
            }
        return s;
    }



    // Painfull key to key name
    std::vector<std::string> GetKeys(std::vector<olc::Key> kv, bool bUppercase = true) {
        std::vector<std::string> s;
        for(auto &k : kv) {
            std::string sTemp = GetKey(k);
            if(sTemp.size() > 0) s.push_back(sTemp);
            if(sTemp.size() > 0 && !bUppercase) s[s.size()-1] = ASCIIToLower(s[s.size()-1].c_str()[0]);
        }
        return s;
    }



    // Painfull key to key name
    std::vector<std::string> GetKeys(std::unordered_set<olc::Key> kv, bool bUppercase = true) {
        std::vector<std::string> s;
        for(auto &k : kv) {
            std::string sTemp = GetKey(k);
            if(sTemp.size() > 0) s.push_back(sTemp);
            if(sTemp.size() > 0 && !bUppercase) s[s.size()-1] = ASCIIToLower(s[s.size()-1].c_str()[0]);
        }
        return s;
    }



    // Gets Caps Lock state crossplatform
    bool GetCapsLockState() {
        #if defined _WIN32
            return ((GetKeyState(VK_CAPITAL) & 0x0001)!=0);
        #elif defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
            X11::Display* disp = X11::XOpenDisplay(NULL);
            // X11::Display *disp = X11::XOpenDisplay(":0");
            X11::XKeyboardState x;
            XGetKeyboardControl(disp, &x);
            X11::XCloseDisplay(disp);
            return (x.led_mask & 1);
        #endif
    }



    // Clipboard for Linux
    #if defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
    char * XPasteType(X11::Atom atom, X11::Display* display, X11::Window window) {
        X11::Atom UTF8 = X11::XInternAtom(display, "UTF8_STRING", True);
        X11::XEvent event;
        int format;
        unsigned long N, size;
        char * data, * s = 0;
        X11::Atom target,
            CLIPBOARD = X11::XInternAtom(display, "CLIPBOARD", 0),
            XSEL_DATA = X11::XInternAtom(display, "XSEL_DATA", 0);
        X11::XConvertSelection(display, CLIPBOARD, atom, XSEL_DATA, window, CurrentTime);
        X11::XSync(display, 0);
        X11::XNextEvent(display, &event);
        
        switch(event.type) {
            case SelectionNotify:
            if(event.xselection.selection != CLIPBOARD) break;
            if(event.xselection.property) {
                XGetWindowProperty(event.xselection.display, event.xselection.requestor,
                    event.xselection.property, 0L,(~0L), 0, AnyPropertyType, &target,
                    &format, &size, &N,(unsigned char**)&data);
                if(target == UTF8 || target == XA_STRING) {
                    s = strndup(data, size);
                    X11::XFree(data);
                }
                XDeleteProperty(event.xselection.display, event.xselection.requestor, event.xselection.property);
            }
        }
        return s;
    }
    #endif



    // Wraper for XPasteType
    std::string XPaste() {
        #if defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
            const char * c = 0;
            X11::Display * display;
            X11::Window window;
            display = X11::XOpenDisplay(0);
            window = X11::XCreateSimpleWindow(display, X11::XDefaultRootWindow(display), 0, 0, 1, 1, 0,0,1);	
            X11::Atom UTF8 = X11::XInternAtom(display, "UTF8_STRING", True);
            if(UTF8 != None) c = XPasteType(UTF8,display,window);
            if(!c) c = XPasteType(UTF8,display,window);
            return std::string(c);
        #endif
        throw "Using Linux function on Windows!";
    }



    // Get Clipboard Content (Crossplatform)
    std::string GetClipboard() {
        #if defined _WIN32
            HANDLE clip;
            if (OpenClipboard(NULL)) {
                clip = GetClipboardData(CF_TEXT);
                CloseClipboard();
            }
            return (char*)clip;
        #elif defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
            return XPaste();
        #endif
    }



    bool RectVsRect(const olc::vf2d Rect1Pos, const olc::vf2d Rect1Size, const olc::vf2d Rect2Pos, const olc::vf2d Rect2Size)
	{
		return (Rect1Pos.x < Rect2Pos.x + Rect2Size.x && Rect1Pos.x + Rect1Size.x > Rect2Pos.x && Rect1Pos.y < Rect2Pos.y + Rect2Size.y && Rect1Pos.y + Rect1Size.y > Rect2Pos.y);
	}



    #pragma endregion Functions




    #pragma region Classes

    

    //  Base class of the XOLCTS GUI system.
    class Renderable 
    {
        public:



            olc::vf2d vLocation = olc::vf2d(0,0);                       //  Vector Storing the location of the object.

            bool bIsChild = false;                                      //  Used to track if an object is a child of another.
            
            bool bActive = true;                                        //  Can be used to create your own enabled / disable system.
            
            int iWidth = 100;                                           //  Width of the object.
            
            int iHeight = 16;                                           //  Height of the object. 
            
            int iType = 0;                                              //  Integer storing the type of object as an ID.
            
            bool bRequireUpdate;                                        //  Indicates the object it needs to be updated (if a custom routine is in the object).

            bool bSelectable = true;                                    //  Tells wether the Renderable can be selected or highlighted (if any interaction is possible). 

            int iHeightChange = 0;                                      // Used by require update to indicate how much to grow

            bool bUpdateWithParent = true;

            virtual void Initialize() {                                 // Routine used to initialize objects 

                std::cout << "Renderable Initialised." << std::endl;    //  Placeholder

                return;

            };



            virtual void UpdateSelf(olc::PixelGameEngine *pge) {        // Routine to update object.

                
                std::cout << "Renderable Updated" << std::endl;         //  Placeholder
                
                return;
            
            };



            virtual void DrawSelf(olc::PixelGameEngine *pge) {          //  Routine to draw object.

                
                std::cout << "Renderable Drawed" << std::endl;          //  Placeholder
                
                return;

            };


        public:



            virtual int GetType() {         //  Returns the type of the object. Used to sort objects from a Renderable array

                return iType;
 
            };



            virtual int GetWidth() {        // Returns the width of the object.

                return iWidth;

            };



            virtual int GetHeight() {       //  Returns the height of the object

                return iHeight;

            };



            virtual void UpdateValue(olc::PixelGameEngine *pge) {       //  Updates the destination Variable

                return;

            };



    };



    // A very basic slider class.
    class Slider : public Renderable 
    {
        public:
            bool bActive;
            // int iWidth = 100;
            // int iHeight = 16;
            int iType = 1;
            int iInitialValue = 0;

            Slider(olc::vf2d VPos, int IInitialValue, float *Variable, std::string Name) {
                vLocation = VPos;
                iSliderMinValue = 0;
                iSliderMaxValue = 100;
                iSliderLenght = 100;
                sSliderText = Name;
                fSliderValue = IInitialValue;
                fDestinationVariable = Variable;
                *fDestinationVariable = IInitialValue;
                iInitialValue = IInitialValue;
            };

            Slider(olc::vf2d VPos, int IMin, int IMax, int ILenght, int IInitialValue, float *Variable) {
                vLocation = VPos;
                iSliderMinValue = IMin;
                iSliderMaxValue = IMax;
                iSliderLenght = ILenght;
                fSliderValue = IInitialValue;
                fDestinationVariable = Variable;
                *fDestinationVariable = IInitialValue;
                iInitialValue = IInitialValue;
            };

            Slider(olc::vf2d VPos, int IMin, int IMax, int ILenght, int IInitialValue, float *Variable, std::string Name) {
                vLocation = VPos;
                iSliderMinValue = IMin;
                iSliderMaxValue = IMax;
                sSliderText = Name;
                iSliderLenght = ILenght;
                fSliderValue = IInitialValue;
                fDestinationVariable = Variable;
                *fDestinationVariable = IInitialValue;
                iInitialValue = IInitialValue;
            };

            Slider(olc::vf2d VPos, int IMin, int IMax, int ILenght, int IInitialValue, std::string Name) {
                vLocation = VPos;
                iSliderMinValue = IMin;
                iSliderMaxValue = IMax;
                sSliderText = Name;
                iSliderLenght = ILenght;
                fSliderValue = IInitialValue;
                iInitialValue = IInitialValue;
            };

            Slider(olc::vf2d VPos, int IMin, int IMax, int ILenght, int IInitialValue) {
                vLocation = VPos;
                iSliderMinValue = IMin;
                iSliderMaxValue = IMax;
                iSliderLenght = ILenght;
                fSliderValue = IInitialValue;
                *fDestinationVariable = IInitialValue;
                iInitialValue = IInitialValue;
            };

            Slider() {
                vLocation = olc::vf2d(0,0);
                iSliderMinValue = 0;
                iSliderMaxValue = 100;
                iSliderLenght = 100;
                iInitialValue = 0;
                // fSliderValue = 0;
                // *fDestinationVariable = fSliderValue;
            };

            virtual void Initialize() {
                *fDestinationVariable = fSliderValue;
                std::cout << "Slider Initialised" << std::endl;
                return;
            };

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                olc::HWButton MouseButttonStateArray[3] {pge->GetMouse(0),pge->GetMouse(1),pge->GetMouse(2)};
                float felapsedtime = pge->GetElapsedTime();
                olc::vf2d MousePos = pge->GetMousePos();
                if(bPressed) {bPressed = false;}
                if(bReleased) {bReleased = false;}
                if(bUpdated) {bUpdated = false;}
                if(XOLCTS::InBounds(vLocation,olc::vf2d(vLocation.x+iWidth+16,vLocation.y+16),MousePos) || bHeld) {
                    bHoovered = true;
                    if(MouseButttonStateArray[1].bPressed) {
                        fSliderValue = iInitialValue;
                    }
                    if(MouseButttonStateArray[0].bPressed) {
                        bPressed = true;
                    }
                    if(MouseButttonStateArray[0].bPressed || bHeld) {
                        bUsesRenderable = true;
                        fSliderValue = Map(MousePos.x-vLocation.x,0,iWidth,iSliderMinValue,iSliderMaxValue);
                        bUpdated = true;
                        bHeld = true;
                    }
                } else {
                    bHoovered = false;
                }
                if(!MouseButttonStateArray[0].bHeld && bHeld) {bReleased = true;bHeld = false;bUpdated = true;}
                if(fSliderValue > iSliderMaxValue) {fSliderValue = iSliderMaxValue;}
                if(fSliderValue < iSliderMinValue) {fSliderValue = iSliderMinValue;}
                *fDestinationVariable = fSliderValue;

                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                std::string sString = sSliderText+" ";
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << fSliderValue;
                sString += std::string("(") + stream.str() + std::string(")");
                pge->SetPixelMode(olc::Pixel::ALPHA);
                pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),olc::Pixel(0,0,0,200));
                pge->FillRect(vLocation,olc::vf2d(Map(fSliderValue,iSliderMinValue,iSliderMaxValue,0,iWidth),iHeight),olc::Pixel(180,164,0));
                pge->DrawString(vLocation+olc::vf2d(4,4),sString,olc::Pixel(255,255,255));
                pge->SetPixelMode(olc::Pixel::NORMAL);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
            virtual void UpdateValue(olc::PixelGameEngine *pge) {
                *fDestinationVariable = fSliderValue;
                return;
            };

        public:

            void SetDestinationVariable(float *Variable) {
                fDestinationVariable = Variable;
                *fDestinationVariable = fSliderValue;
            }

        public:
            float fSliderValue = 0;
            float* fDestinationVariable = &fSliderValue;
            bool bHeld = false;
            bool bPressed = false;
            bool bReleased = false;
            bool bHoovered = false;
            bool bUpdated = false;
            std::string sSliderText = "Slider";
            int InputCooldown = 0;
            int iSliderLenght;
            int iSliderMaxValue;
            int iSliderMinValue;
    };



    // Simple two states button.
    class Button : public Renderable
    {
        public:
            bool bActive;
            int iType = 3;

            Button(olc::vf2d VPos) {
                vLocation = VPos;
                *bDestinationVariable = bButtonValue;
            };            
            Button(olc::vf2d VPos, bool KeepValue, std::string name, bool *Variable) {
                bDestinationVariable = Variable;
                vLocation = VPos;
                *bDestinationVariable = bButtonValue;
                bKeepValue = KeepValue;
                sButtonText = name;
            };
            
            Button(olc::vf2d VPos, bool KeepValue, std::string name) {
                vLocation = VPos;
                *bDestinationVariable = bButtonValue;
                bKeepValue = KeepValue;
                sButtonText = name;
            };

            Button() {
                vLocation = olc::vf2d(0,0);
                *bDestinationVariable = bButtonValue;
            };

            virtual void Initialize() {
                std::cout << "Renderable Initialised." << std::endl;
                return;
            };

            virtual bool CustomUpdate() 
            {
                return false;;
            };

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                olc::HWButton MouseButttonStateArray[3] {pge->GetMouse(0),pge->GetMouse(1),pge->GetMouse(2)};
                olc::vf2d MousePos = pge->GetMousePos();
                bPressed = false;
                bReleased = false;
                bReleasedOnButton = false;
                if(!bKeepValue) {bButtonValue = false;}
                if(XOLCTS::InBounds(vLocation+olc::vf2d((iWidth-16)*!(!bKeepValue && !bShowRect),0),olc::vf2d(vLocation.x+iWidth,vLocation.y+16),MousePos) || bHeld) {
                    bHoovered = true;
                    if(MouseButttonStateArray[0].bPressed) {
                        bUsesRenderable = true;
                        bPressed = true;
                    }
                    if(MouseButttonStateArray[0].bPressed || (bHeld && !bKeepValue)) {
                        bButtonValue = !bButtonValue;
                        bHeld = true;
                    }
                } else {
                    bHoovered = false;
                }
                if(!MouseButttonStateArray[0].bHeld && bHeld) {
                    bReleased = true;
                    bHeld = false;
                    if(XOLCTS::InBounds(vLocation+olc::vf2d((iWidth-16)*!(!bKeepValue && !bShowRect),0),olc::vf2d(vLocation.x+iWidth,vLocation.y+16),MousePos)) bReleasedOnButton = true;
                }
                *bDestinationVariable = bButtonValue;
                CustomUpdate();
                return;
            };


        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
            virtual void UpdateValue(olc::PixelGameEngine *pge) {
                *bDestinationVariable = bButtonValue;
                return;
            };

        public:

            void SetDestinationVariable(bool *Variable) {
                bDestinationVariable = Variable;
                *bDestinationVariable = bButtonValue;
            }

        public:
            bool bButtonValue = false;
            bool bKeepValue = true;
            bool* bDestinationVariable = &bButtonValue;
            bool bHeld = false;
            bool bPressed = false;
            bool bReleased = false;
            bool bHoovered = false;
            bool bShowRect = false;
            bool bReleasedOnButton = false;
            olc::Pixel pTextColor = olc::WHITE;
            std::string sButtonText = "Button";
    };



    // Simple two states button.
    class ButtonArray : public Renderable
    {
        public:
            bool bActive;
            int iWidth = 140;
            int iHeight = 16;
            int iType = 5;

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {

                if(!bIsChild) {     // Allows moving menu
                    if(pge->GetMouse(0).bReleased) {bHoldButtonArray = false;vOffset = olc::vf2d(0,0);}
                    if(XOLCTS::InBounds(olc::vf2d(vLocation.x,vLocation.y),olc::vf2d(vLocation.x+iWidth-16,vLocation.y+16), pge->GetMousePos()) && pge->GetMouse(0).bPressed) {
                        vOffset = pge->GetMousePos()-vLocation;
                        bHoldButtonArray = true;
                        return;
                    }
                    if(bHoldButtonArray) {
                        bUsesRenderable = true;
                        vLocation += pge->GetMousePos()-vLocation-vOffset;
                        bRequireUpdate = true;
                    }
                }

                if(bRequireUpdate) {
                    iButtonArrayHeight = 16;
                    iHeight = 32-(16*(1-bExpanded));
                    if(bExpanded) {
                        for(auto &i : vButtons) {
                            i->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iButtonArrayHeight+1);
                            i->bRequireUpdate = true;
                            i->UpdateSelf(pge);
                            iButtonArrayHeight += i->GetHeight()+1;
                            iHeight += i->GetHeight()+1;
                        } 
                    } 
                    bRequireUpdate = false;
                    return;
                }



                if(XOLCTS::InBounds(olc::vf2d(vLocation.x+(iWidth-16),vLocation.y),olc::vf2d(vLocation.x+iWidth,vLocation.y+16), pge->GetMousePos()) && pge->GetMouse(0).bPressed) {
                    bExpanded = !bExpanded;
                    bUsesRenderable = true;
                    iButtonArrayHeight = 16;
                    iHeight = 32-(16*(1-bExpanded));
                    if(bExpanded) {
                        for(auto &i : vButtons) {
                            i->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iButtonArrayHeight+1);
                            i->UpdateSelf(pge);
                            iButtonArrayHeight += i->GetHeight()+1;
                            iHeight += i->GetHeight()+1;
                        } 
                    } 
                    return;
                }
                if(!bExpanded) return;



                iButtonArrayHeight = 16;
                iHeight = 32-(16*(1-bExpanded));

                int iIndex = 0;
                bool bAnyButtonPressed = false;
 
                for(auto &i : vButtons) {
                    i->UpdateSelf(pge);
                    if(i->bButtonValue) bAnyButtonPressed = true;
                    if(i->bButtonValue && i->bPressed) {
                        for(auto &j : vButtons) {
                            j->bButtonValue = false;
                        }
                        i->bButtonValue = true;
                        iIndexOfActiveButton = iIndex;
                    }
                    iButtonArrayHeight += i->GetHeight()+1;
                    iHeight += i->GetHeight()+1;
                    iIndex++;
                }

                if(!bAnyButtonPressed) iIndexOfActiveButton = -1;
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                pge->SetPixelMode(olc::Pixel::ALPHA);
                if(!bIsChild) {
                    pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight*bExpanded+(16*(1-bExpanded))),olc::Pixel(128,128,128,200));
                }
                pge->DrawRect(vLocation,olc::vf2d(0,iHeight*bExpanded+(16*(1-bExpanded))),olc::Pixel(200,200,200,100));
                pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight*bExpanded+(16*(1-bExpanded))),olc::Pixel(128,128,128,64));
                pge->FillRect(vLocation,olc::vf2d(iWidth,16),olc::Pixel(50,50,50,128));
                // pge->FillRect(vLocation+olc::vf2d(iWidth-16,0),olc::vf2d(16,16),olc::Pixel(255,255,255,128));
                pge->DrawRect(vLocation+olc::vf2d(iWidth-16,0),olc::vf2d(15,15),olc::Pixel(255,255,255,128));
                pge->DrawString(vLocation+olc::vf2d(4,4),sButtonArrayName);
                pge->SetPixelMode(olc::Pixel::NORMAL);
                if(!bExpanded) return;
                for(auto &i : vButtons) {
                    i->DrawSelf(pge);
                }
                return;
            };

            void SetValue(int value) {
                if(value < vButtons.size()) {
                    iIndexOfActiveButton = value;
                    for(auto &i : vButtons) {
                        i->bButtonValue = false;
                    }
                    vButtons[value]->bButtonValue = true;
                }
                return;
            }

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };


        public:
            std::vector<XOLCTS::Button*> vButtons;
            bool bExpanded = false;
            std::string sButtonArrayName = "Button Array";
            int iButtonArrayHeight = 16;
            // bool bOneChoice = true;
            int iIndexOfActiveButton = -1;

            void AddElement(XOLCTS::Button* e) {
                e->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iButtonArrayHeight+1);
                e->bIsChild = true;
                *e->bDestinationVariable = e->bButtonValue;
                vButtons.push_back(e);
            };

        private:
            bool bHoldButtonArray = false;
            olc::vf2d vOffset;
    };



    // Basic menu class. Used to list other Renderables.
    class Menu : public Renderable 
    {        
        public:
            bool bActive;
            int iWidth = 132;
            int iHeight = 16;
            int iType = 2;

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                if(!bIsChild) {
                    if(pge->GetMouse(0).bReleased) {bHoldMenu = false;vOffset = olc::vf2d(0,0);}
                    if(bMovable && XOLCTS::InBounds(olc::vf2d(vLocation.x,vLocation.y),olc::vf2d(vLocation.x+iWidth-16,vLocation.y+16), pge->GetMousePos()) && pge->GetMouse(0).bPressed) {
                        vOffset = pge->GetMousePos()-vLocation;
                        bHoldMenu = true;
                        return;
                    }
                    if(bHoldMenu) {
                        vLocation += pge->GetMousePos()-vLocation-vOffset;
                        bRequireUpdate = true;
                    }
                }

                if(bRequireUpdate) { 
                    iMenuHeight = 16;
                    iHeight = 32-(16*(1-bExpanded));
                    if(bExpanded) {
                        for(auto &i : vElements) {
                            i->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iMenuHeight+1);
                            i->bRequireUpdate = true;
                            i->UpdateSelf(pge);
                            iMenuHeight += i->GetHeight()+1;
                            iHeight += i->GetHeight()+1;
                        } 
                    } 
                    bRequireUpdate = false;
                    return;
                }



                if(XOLCTS::InBounds(olc::vf2d(vLocation.x+(iWidth-16),vLocation.y),olc::vf2d(vLocation.x+iWidth,vLocation.y+16), pge->GetMousePos()) && pge->GetMouse(0).bPressed) {
                    bExpanded = !bExpanded;
                    bUsesRenderable = true;
                    // iMenuHeight = 16;
                    // iHeight = 32-(16*(1-bExpanded));
                    // if(bExpanded) {
                    //     for(auto &i : vElements) {
                    //         i->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iMenuHeight+1);
                    //         i->UpdateSelf(pge);
                    //         iMenuHeight += i->GetHeight()+1;
                    //         iHeight += i->GetHeight()+1;
                    //     } 
                    // } 
                    // return;
                }
                if(!bExpanded) return;

                iMenuHeight = 16;
                iWidth = 132;
                iHeight = 32-(16*(1-bExpanded));

                int record = 0;

                for(auto &i : vElements) {
                    i->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iMenuHeight+1);
                    // if(i->GetType() == 2 || i->GetType() == 5) i->iWidth += 8;
                    i->UpdateSelf(pge);
                    if(i->GetWidth() > record) {
                        record = i->GetWidth();
                    }
                    // if(i->GetType() == 2 || i->GetType() == 5) iWidth -= i->iWidth-iWidth;
                    // iWidth += i->GetWidth();
                    iMenuHeight += i->GetHeight()+1;
                    iHeight += i->GetHeight()+1;
                }
                iWidth = record+16;

                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                pge->SetPixelMode(olc::Pixel::ALPHA);
                if(!bIsChild) {
                    pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight*bExpanded+(16*(1-bExpanded))),olc::Pixel(100,100,100,220));
                }
                pge->DrawRect(vLocation,olc::vf2d(0,iHeight*bExpanded+(16*(1-bExpanded))),olc::Pixel(200,200,200,100));
                pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight*bExpanded+(16*(1-bExpanded))),olc::Pixel(128,128,128,64));
                pge->FillRect(vLocation,olc::vf2d(iWidth,16),olc::Pixel(50,50,50,128));
                pge->DrawRect(vLocation+olc::vf2d(iWidth-16,0),olc::vf2d(15,15),olc::Pixel(255,255,255,128));
                pge->DrawString(vLocation+olc::vf2d(4,4),sMenuName);
                pge->SetPixelMode(olc::Pixel::NORMAL);
                if(!bExpanded) return;
                for(auto &i : vElements) {
                    i->DrawSelf(pge);
                }
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
            
            virtual void UpdateValue(olc::PixelGameEngine *pge) {
                for(auto &i : vElements) {
                    i->UpdateValue(pge);
                    i->UpdateSelf(pge);
                }
                return;
            }


        public:
            std::vector<XOLCTS::Renderable*> vElements;
            bool bExpanded = false;
            std::string sMenuName = "Menu";
            int iMenuHeight = 16;
            bool bMovable = true;

            void ResetMenu() {
                iMenuHeight = 16;
                iWidth = 132;
                for(auto &i : vElements) {
                    delete i;
                }
                vElements.clear();
            }

            int AddElement(XOLCTS::Renderable* e) {
                e->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iMenuHeight+1);
                e->bIsChild = true;
                if(e->GetType() == 2) iWidth += 8;
                if(e->GetType() == 5) iWidth += 8;
                vElements.push_back(e);
                
                iWidth = 132;

                int record = 0;

                    if(e->GetWidth() > record) {
                        record = e->GetWidth();
                    }

                iWidth = record+8;
                return vElements.size()-1;
            };

        private:
            bool bHoldMenu = false;
            olc::vf2d vOffset;
    };



    // Simple separator for menus
    class Separator : public Renderable  
    {
        public:
            // int iWidth = 100;
            // int iHeight = 8;
            int iType = 8;
            // bool bSelectable = false;

            Separator(olc::vf2d vPos, int Width) {
                iWidth = Width;
                vLocation = vPos;
                bSelectable = false;
                iHeight = 9;
            }

            Separator() {
                vLocation = olc::vf2d(0,0);
                bSelectable = false;
                iHeight = 9;
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                pge->DrawLine(vLocation+olc::vf2d(0,iHeight/2), vLocation+olc::vf2d(iWidth,iHeight/2), olc::Pixel(128,128,128));
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
    };



    // Basic menu class. Used to list other Renderables.
    class SubMenu : public Renderable 
    {        
        public:
            bool bActive;
            // int iHeight = 16;
            int iType = 9;


            SubMenu() {
                return;
            };

            SubMenu(olc::vf2d Location, int Width, int MenuWidth, std::string Text, int MenuID) {
                sMenuName = Text;
                iMenuID = MenuID;
                iWidth = Width;
                iMenuWidth = MenuWidth; 
                vLocation = Location;
            };

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {

                // Variables Update

                if(iCloseTimer > 0) iCloseTimer--;
                if(iCloseTimer < 1 ) iCloseTimer = 0;
                if(iSubMenuOpened != iMenuID) iCloseTimer = 0;


                // Basicly Getting the hovered menu element and store it to draw highlight later.
                if(iSubMenuOpened == iMenuID) {
                    bUsesRenderable = true;
                    int index = 0;
                    iHoveredElement = -1;
                    for(auto &i : vElements) {
                        i->UpdateSelf(pge);
                        if(i->bSelectable && InBounds(i->vLocation-olc::vf2d(8,0),i->vLocation+olc::vf2d(iMenuWidth,i->iHeight),pge->GetMousePos())) iHoveredElement = index;
                        index++;
                    }
                }
                
                // Inputs Handling.
                if(InBounds(vLocation,vLocation+olc::vf2d(iWidth+8,16),pge->GetMousePos())) {iCloseTimer = 60; Open(); return;} // If a menu is opened open the one hovered else if none is opened open the one clicked.
                if(iCloseTimer > 0 && InBounds(vLocation+olc::vf2d(iWidth+8,-2),vLocation+olc::vf2d(iWidth+iMenuWidth,iMenuHeight+7),pge->GetMousePos())) {iCloseTimer = 60; Open(); return;} // If a menu is opened open the one hovered else if none is opened open the one clicked.
                if(iSubMenuOpened == iMenuID && iCloseTimer == 0) iSubMenuOpened = -1; // Closes menu if clicked on empty space.
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                pge->SetPixelMode(olc::Pixel::ALPHA);
                pge->DrawString(vLocation+olc::vf2d(0,4),sMenuName,olc::Pixel(200,200,200)); // Draws menu text (name).
                pge->FillTriangle(olc::vf2d(vLocation.x+iWidth-8,vLocation.y+6),olc::vf2d(vLocation.x+iWidth-2,vLocation.y+(iHeight/2)),olc::vf2d(vLocation.x+iWidth-8,vLocation.y+iHeight-6), olc::Pixel(160,160,160));
                if(iCloseTimer > 0 || iSubMenuOpened == iMenuID) { // Stops drawing if menu isn't opened (prevents from drawing the menu elements).
                    pge->FillRect(vLocation+olc::vf2d(iWidth+8,-2),olc::vf2d(iMenuWidth+8,iMenuHeight+7),pBackgroundColor); // Draws menu background box.
                    pge->DrawRect(vLocation+olc::vf2d(iWidth+8,-2),olc::vf2d(iMenuWidth+8,iMenuHeight+7),olc::BLACK); // Draws menu box outline.
                    if(iHoveredElement != -1) pge->FillRect(vElements[iHoveredElement]->vLocation-olc::vf2d(7,0),olc::vf2d(iMenuWidth+6,vElements[iHoveredElement]->iHeight),olc::Pixel(50, 50, 50)); // Draws element highlight if necessary.

                    pge->SetPixelMode(olc::Pixel::NORMAL);
                    for(auto &i : vElements) {
                        i->DrawSelf(pge);
                    }
                    return;
                } 
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };


        public:
            std::vector<XOLCTS::Renderable*> vElements;
            std::string sMenuName = "Menu";
            olc::Pixel pColor = olc::Pixel(37,37,37);
            int iMenuWidth = 0;
            int iMenuHeight = 2;
            int iCloseTimer = 0; 
            int iHoveredElement = -1; 
            int iMenuID;

            void Open() {
                if(iSubMenuOpened == iMenuID) return;
                iSubMenuOpened = iMenuID;
                iCloseTimer = 60;
            }

            void AddElement(XOLCTS::Renderable* e) {
                e->vLocation = olc::vf2d(vLocation.x+iWidth+16,vLocation.y+iMenuHeight+1);
                e->bIsChild = true;
                e->iWidth = iMenuWidth-16;
                iMenuHeight += e->iHeight;
                vElements.push_back(e);
            };

            void AddSeparator() {
                AddElement(new Separator(olc::vf2d(0,0),iMenuWidth-16));
            }
    };



    // Basic menu class. Used to list other Renderables.
    class TopMenu : public Renderable 
    {        
        public:
            // bool bActive;
            int iWidth = 0;
            int iHeight = 16;
            int iType = 7;


            TopMenu() {
                return;
            };

            TopMenu(olc::vf2d Location, int Width, int MenuWidth, std::string Text, int MenuID) {
                sMenuName = Text;
                iMenuID = MenuID;
                iWidth = Width;
                iMenuWidth = MenuWidth; 
                vLocation = Location;
            };

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {

                // Variables Update
                if(!bActive) return;
                if(iMenuOpened == iMenuID) { // Variable for menu white highlight Width.
                    bUsesRenderable = true;
                    if(iHighlighterWidth < iWidth/2) iHighlighterWidth += iWidth/16; 
                    if(iHighlighterWidth > iWidth/2) iHighlighterWidth = iWidth/2;
                }
                if(iMenuOpened != iMenuID) { // Variable for menu white highlight Opacity.
                    if(iHighlighterOpacity > 1) iHighlighterOpacity -= 255/16;
                    if(iHighlighterWidth < 1) iHighlighterWidth = 0;
                } else {
                    iHighlighterOpacity = 255;
                }


                // Basicly Getting the hovered menu element and store it to draw highlight later.
                if(iMenuOpened == iMenuID) {
                    int index = 0;
                    iHoveredElement = -1;
                    for(auto &i : vElements) {
                        if(!i->bUpdateWithParent || !i->bActive) continue;
                        i->UpdateSelf(pge);
                        if(i->bRequireUpdate) {iMenuHeight += i->iHeightChange; i->bRequireUpdate = false; i->iHeightChange = 0;};
                        if(i->bSelectable && InBounds(i->vLocation-olc::vf2d(8,0),i->vLocation+olc::vf2d(iMenuWidth,i->iHeight),pge->GetMousePos())) iHoveredElement = index;
                        index++;
                    }
                }
                
                // Inputs Handling.
                if(iMenuOpened != -1 || pge->GetMouse(0).bPressed) {
                    if(InBounds(vLocation,vLocation+olc::vf2d(iWidth,iHeight),pge->GetMousePos())) { Open(); return;} // If a menu is opened open the one hovered else if none is opened open the one clicked.
                    if(iMenuOpened == iMenuID && InBounds(vLocation+olc::vf2d(0,iHeight),vLocation+olc::vf2d(iMenuWidth,iMenuHeight+iHeight),pge->GetMousePos())) {iMenuOpened = iMenuID; return;} // Just detects if the user clicked in the menu box to prevent menu closing.
                    if(iMenuOpened == iMenuID && iSubMenuOpened == -1 && pge->GetMouse(0).bPressed)iMenuOpened = -1; // Closes menu if clicked on empty space.
                    if(pge->GetKey(olc::ESCAPE).bPressed) iMenuOpened = -1; // Closes menu if escape is pressed.
                }

                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                if(!bActive) return;
                pge->SetPixelMode(olc::Pixel::ALPHA);
                pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),pColor); // Drawing the menu text box.
                if(iHighlighterOpacity > 0) pge->FillRect(vLocation+olc::vf2d(iWidth/2-iHighlighterWidth,0),olc::vf2d(iHighlighterWidth*2,iHeight),olc::Pixel(64,64,64,iHighlighterOpacity)); // Draws menu text box highlight
                pge->DrawString(vLocation+olc::vf2d(0,4),sMenuName,olc::Pixel(200,200,200)); // Draws menu text (name).
                if(iMenuOpened != iMenuID) return; // Stops drawing if menu isn't opened (prevents from drawing the menu elements).
                pge->FillRect(vLocation+olc::vf2d(0,iHeight),olc::vf2d(iMenuWidth,iMenuHeight+3),pBackgroundColor); // Draws menu background box.
                pge->DrawRect(vLocation+olc::vf2d(0,iHeight),olc::vf2d(iMenuWidth,iMenuHeight+3),olc::BLACK); // Draws menu box outline.
                // if(iHoveredElement != -1) pge->FillRect(vElements[iHoveredElement]->vLocation-olc::vf2d(8,0),olc::vf2d(iMenuWidth,vElements[iHoveredElement]->iHeight),olc::Pixel(39, 116, 252)); // Draws element highlight if necessary.
                if(iHoveredElement != -1) pge->FillRect(vElements[iHoveredElement]->vLocation-olc::vf2d(7,0),olc::vf2d(iMenuWidth-2,vElements[iHoveredElement]->iHeight),olc::Pixel(50, 50, 50)); // Draws element highlight if necessary.

                pge->SetPixelMode(olc::Pixel::NORMAL);
                for(auto &i : vElements) {
                    i->DrawSelf(pge);
                    pge->SetPixelMode(olc::Pixel::ALPHA);
                    if(!i->bActive) pge->FillRect(i->vLocation,olc::vf2d(i->iWidth,i->iHeight),olc::Pixel(pBackgroundColor.r/2,pBackgroundColor.g/2,pBackgroundColor.b/2,223));
                    pge->SetPixelMode(olc::Pixel::NORMAL);
                }
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };


        public:
            std::vector<XOLCTS::Renderable*> vElements;
            std::string sMenuName = "Menu";
            olc::Pixel pColor = olc::Pixel(37,37,37);
            int iMenuWidth = 0;
            int iMenuHeight = 2;
            int iSubMenuTimer = 0;
            int iHighlighterWidth = 0; 
            int iHighlighterOpacity = 0;
            int iHoveredElement = -1; 
            int iMenuID;
            // olc::Pixel pBackgroundColor = olc::Pixel(35,35,35);
            // olc::Pixel pBackgroundColor = olc::Pixel(38,37,37);

            void Open() {
                if(iMenuOpened == iMenuID) return;
                iMenuOpened = iMenuID;
                // std::cout << std::to_string(iMenuID)+" Opened." << std::endl;
                iHighlighterWidth = 0; 
                iSubMenuTimer = 0;
                iHighlighterOpacity = 255;
            }

            void AddElement(XOLCTS::Renderable* e) {
                e->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iHeight+iMenuHeight+1);
                e->bIsChild = true;
                e->iWidth = iMenuWidth-16;
                iMenuHeight += e->iHeight+1;
                vElements.push_back(e);
            };

            void AddSeparator() {
                AddElement(new Separator(olc::vf2d(0,0),iMenuWidth-16));
            }

            SubMenu* AddSubmenu(olc::vf2d Location, int Width, int MenuWidth, std::string Text, int MenuID) {
                SubMenu* i = new SubMenu(Location, Width, MenuWidth, Text, MenuID);
                AddElement(i);
                return  i;
            }

            template <typename T>
            std::vector<T> GetAllOfClass() {
                std::vector<T> v;
                for(auto &i : vElements) {
                    if(dynamic_cast<T>(i) != nullptr) v.push_back(dynamic_cast<T>(i));
                }
                return v;
            }

            template <typename T>
            std::vector<int> GetAllIndexOfClass() {
                int c = 0;
                std::vector<int> v;
                for(auto &i : vElements) {
                    if(dynamic_cast<T>(i) != nullptr) v.push_back(c);
                    c++;
                }
                return v;
            }

            int GetIndexOfItem(Renderable* r) {
                int c = 0;
                for(auto &i : vElements) {
                    if(r == i) return c;
                    c++;
                }
                return -1;
            }

            void ResetMenu() {
                for(auto &i : vElements) {
                    delete i;
                }
                vElements.clear();
                
                iMenuHeight = 2;
                // int iHeight = 16;
            }
    };



    // Display text object.
    class Text : public Renderable  
    {
        public:
            // int iWidth = 100;
            // int iHeight = 16;
            int iType = 4;

            Text(olc::vf2d vPos) {
                vLocation = vPos;
            }

            Text(olc::vf2d vPos, std::string s) {
                sText = s;
                vLocation = vPos;
            }

            Text() {
                vLocation = olc::vf2d(0,0);
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                pge->DrawString(vLocation+olc::vf2d(0,(iHeight-8)/2),sText,pColor);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        
        public:
            std::string sText = "Text";
            olc::Pixel pColor = olc::Pixel(128,128,128);
    };



    // Simple Text Input
    class TextInput : public Renderable  
    {
        public:
            // int iWidth = 100;
            // int iHeight = 16;
            int iType = 10;

            TextInput(olc::vf2d vPos) {
                vLocation = vPos;
                bEditable = true;
            }

            TextInput(olc::vf2d vPos, std::string sText) {
                sInputText = sText;
                vLocation = vPos;
                bEditable = true;
            }

            TextInput() {
                vLocation = olc::vf2d(0,0);
                bEditable = true;
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                iMaxCharsIndex = floor(iWidth/8);
                if(!bEditable) return;
                if(pge->GetMouse(0).bPressed) {
                    bSelected = InBounds(vLocation,vLocation+olc::vf2d(iWidth,iHeight),pge->GetMousePos());
                    iCursorIndex = (int)((pge->GetMousePos().x-vLocation.x)/8)*bSelected;
                }
                if(bSelected) {
                    bUsesRenderable = true;
                    iCursorIndex = iCursorIndex+(1*pge->GetKey(olc::RIGHT).bPressed || (pge->GetKey(olc::RIGHT).bHeld && iKeyTimeout == 0));
                    iCursorIndex = iCursorIndex-(1*pge->GetKey(olc::LEFT).bPressed || (pge->GetKey(olc::LEFT).bHeld && iKeyTimeout == 0));
                    if(pge->GetKey(olc::LEFT).bHeld && iKeyTimeout == 0) iKeyTimeout = 10;
                    if(pge->GetKey(olc::RIGHT).bHeld && iKeyTimeout == 0) iKeyTimeout = 10;

                    if(iCursorIndex < 0) iCursorIndex = 0;
                    if(iCursorIndex > sInputText.size()) iCursorIndex = sInputText.size();
                    if(iCursorIndex-iDisplayOffset > iMaxCharsIndex) iDisplayOffset++;
                    if(iCursorIndex-iDisplayOffset < 1) iDisplayOffset = (iDisplayOffset > iMaxCharsIndex-1) ? iDisplayOffset-iMaxCharsIndex : 0;
                    std::vector<olc::Key> vKeysPressed;
                    for(std::map<size_t, uint8_t>::iterator it = olc::mapKeys.begin(); it != olc::mapKeys.end(); ++it) {
                        if(pge->GetKey(olc::Key(it->second)).bPressed && std::count(vKeysPressed.begin(),vKeysPressed.end(),olc::Key(it->second)) == 0) {vKeysPressed.push_back(olc::Key(it->second));}
                    }

                    std::vector<std::string> sUserInput = GetKeys(vKeysPressed);
                    for(auto &i : sUserInput) {
                        sInputText.insert(iCursorIndex,i);
                        iCursorIndex++; 
                    }
                    iKeyTimeout -= (iKeyTimeout < 1) ? 0 : 1*(pge->GetElapsedTime()*2);
                    if(((pge->GetKey(olc::BACK).bHeld && iKeyTimeout == 0) || pge->GetKey(olc::BACK).bPressed) && iCursorIndex > 0 ) {sInputText.erase(sInputText.begin()+iCursorIndex-1); iCursorIndex--; iKeyTimeout = iKeyTimeoutDuration*2*(pge->GetKey(olc::BACK).bPressed+1);} 
                    if(((pge->GetKey(olc::DEL).bHeld && iKeyTimeout == 0) || pge->GetKey(olc::DEL).bPressed) && iCursorIndex < sInputText.size() && iKeyTimeout == 0) {sInputText.erase(sInputText.begin()+iCursorIndex); iKeyTimeout = iKeyTimeoutDuration*2*(pge->GetKey(olc::DEL).bPressed+1);} 
                    bReturned = pge->GetKey(olc::ENTER).bPressed; 
                    // if(bAdded[olc::RETURN] && iCursorIndex > 0) {sInputText.erase(iCursorIndex); iCursorIndex--; } 
                }
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                std::string sDrawnString = sInputText;
                if(sInputText.size() > iMaxCharsIndex) {
                    sDrawnString = "";
                    for(int i = iDisplayOffset*bSelected; (i < (iDisplayOffset*bSelected)+iMaxCharsIndex && i < sInputText.size()); i++) {
                        sDrawnString += sInputText[i];
                    }
                } else iDisplayOffset = 0;
                
                pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),pBackgroundColor*(!bEditable+1));
                pge->DrawRect(vLocation,olc::vf2d(iWidth,iHeight),olc::Pixel(0,0,0));
                if(sInputText.size() > 0) {
                    pge->DrawString(vLocation+olc::vf2d(4,(iHeight-8)/2),sDrawnString,pColor);
                } else {
                    sHintText = ""; 
                    for(int i = 0; i < iMaxCharsIndex; i++) {
                        sHintText += "-";
                    }
                    pge->DrawString(vLocation+olc::vf2d(4,(iHeight-8)/2),sHintText,pColor);
                }
                if(bSelected) pge->FillRect(vLocation+olc::vf2d(4+(8*(iCursorIndex-iDisplayOffset)),2),olc::vf2d(1,iHeight-2),olc::WHITE);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        
        public:
            std::string sHintText = "";
            std::string sInputText = "";
            bool bSelected = false;
            int iCursorIndex = 0;
            int iDisplayOffset = 0;
            int iMaxCharsIndex = 0;
            int iKeyTimeout = 0;
            int iKeyTimeoutDuration = 3;
            bool bEditable;
            olc::Pixel pColor = olc::Pixel(128,128,128);
            bool bReturned = false;
    };



    // Simple Number Input
    class NumberInput : public Renderable  
    {
        public:
            // int iWidth = 100;
            // int iHeight = 16;
            int iType = 11;

            NumberInput(olc::vf2d vPos) {
                vLocation = vPos;
            }

            NumberInput(olc::vf2d vPos, std::string sText) {
                sInputText = sText;
                vLocation = vPos;
            }

            NumberInput() {
                vLocation = olc::vf2d(0,0);
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                iMaxCharsIndex = floor(iWidth/8);
                if(pge->GetMouse(0).bPressed) {
                    bSelected = InBounds(vLocation,vLocation+olc::vf2d(iWidth,iHeight),pge->GetMousePos());
                    iCursorIndex = (int)((pge->GetMousePos().x-vLocation.x)/8)*bSelected;
                }
                if(bSelected) {
                    bUsesRenderable = true;
                    iCursorIndex = iCursorIndex+(1*pge->GetKey(olc::RIGHT).bPressed || (pge->GetKey(olc::RIGHT).bHeld && iKeyTimeout == 0));
                    iCursorIndex = iCursorIndex-(1*pge->GetKey(olc::LEFT).bPressed || (pge->GetKey(olc::LEFT).bHeld && iKeyTimeout == 0));
                    if(pge->GetKey(olc::LEFT).bHeld && iKeyTimeout == 0) iKeyTimeout = 10;
                    if(pge->GetKey(olc::RIGHT).bHeld && iKeyTimeout == 0) iKeyTimeout = 10;

                    if(iCursorIndex < 0) iCursorIndex = 0;
                    if(iCursorIndex > sInputText.size()) iCursorIndex = sInputText.size();
                    if(iCursorIndex-iDisplayOffset > iMaxCharsIndex) iDisplayOffset++;
                    if(iCursorIndex-iDisplayOffset < 1) iDisplayOffset = (iDisplayOffset > iMaxCharsIndex-1) ? iDisplayOffset-iMaxCharsIndex : 0;
                    std::vector<olc::Key> vKeysPressed;
                    for(std::map<size_t, uint8_t>::iterator it = olc::mapKeys.begin(); it != olc::mapKeys.end(); ++it) {
                        if(pge->GetKey(olc::Key(it->second)).bPressed && std::count(vKeysPressed.begin(),vKeysPressed.end(),olc::Key(it->second)) == 0) {vKeysPressed.push_back(olc::Key(it->second));}
                    }
                    std::vector<std::string> sUserInput = GetKeys(vKeysPressed);
                    for(auto &i : sUserInput) {
                        if(std::isdigit(i[0]) || i == ".") {sInputText.insert(iCursorIndex,i); iCursorIndex++;} 
                    }
                    iKeyTimeout -= (iKeyTimeout < 1) ? 0 : 1*(pge->GetElapsedTime()*2);
                    if(((pge->GetKey(olc::BACK).bHeld && iKeyTimeout == 0) || pge->GetKey(olc::BACK).bPressed) && iCursorIndex > 0 ) {sInputText.erase(sInputText.begin()+iCursorIndex-1); iCursorIndex--; iKeyTimeout = iKeyTimeoutDuration*2*(pge->GetKey(olc::BACK).bPressed+1);} 
                    if(((pge->GetKey(olc::DEL).bHeld && iKeyTimeout == 0) || pge->GetKey(olc::DEL).bPressed) && iCursorIndex < sInputText.size() && iKeyTimeout == 0) {sInputText.erase(sInputText.begin()+iCursorIndex); iKeyTimeout = iKeyTimeoutDuration*2*(pge->GetKey(olc::DEL).bPressed+1);} 
                    bReturned = pge->GetKey(olc::ENTER).bPressed; 
                    // if(bAdded[olc::RETURN] && iCursorIndex > 0) {sInputText.erase(iCursorIndex); iCursorIndex--; } 
                }
                if(sInputText.size()>0) fValue = std::stof(sInputText);
                else fValue = 0;
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                std::string sDrawnString = sInputText;
                if(sInputText.size() > iMaxCharsIndex) {
                    sDrawnString = "";
                    for(int i = iDisplayOffset*bSelected; (i < (iDisplayOffset*bSelected)+iMaxCharsIndex && i < sInputText.size()); i++) {
                        sDrawnString += sInputText[i];
                    }
                } else iDisplayOffset = 0;
                
                pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),pBackgroundColor);
                pge->DrawRect(vLocation,olc::vf2d(iWidth,iHeight),olc::Pixel(0,0,0));
                if(sInputText.size() > 0) {
                    pge->DrawString(vLocation+olc::vf2d(4,(iHeight-8)/2),sDrawnString,pColor);
                } else {
                    sHintText = ""; 
                    for(int i = 0; i < iMaxCharsIndex; i++) {
                        sHintText += "-";
                    }
                    pge->DrawString(vLocation+olc::vf2d(4,(iHeight-8)/2),sHintText,pColor);
                }
                if(bSelected) pge->FillRect(vLocation+olc::vf2d(4+(8*(iCursorIndex-iDisplayOffset)),2),olc::vf2d(1,iHeight-2),olc::WHITE);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        
        public:
            std::string sHintText = "";
            std::string sInputText = "";
            float fValue = 0;
            bool bSelected = false;
            int iCursorIndex = 0;
            int iDisplayOffset = 0;
            int iMaxCharsIndex = 0;
            int iKeyTimeout = 0;
            int iKeyTimeoutDuration = 3;
            olc::Pixel pColor = olc::Pixel(128,128,128);
            bool bReturned = false;
    };



    // Basicly allows to create a line that contains multiple renderables
    class CustomLine : public Renderable 
    {
           public:
            // int iWidth = 100;
            // int iHeight = 16;
            int iType = 12;

            CustomLine(olc::vf2d vPos) {
                vLocation = vPos;
            }

            CustomLine() {
                vLocation = olc::vf2d(0,0);
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                // if(bRequireUpdate) {
                //     bRequireUpdate = false;
                //     for(auto &i : vElements) {
                //         i->vLocation = olc::vf2d(vLocation.x+i->vLocation.x,vLocation.y+i->vLocation.y+1);
                //     }
                // }
                for(auto &i : vElements) {
                    i->vLocation.y = vLocation.y+1;
                    i->UpdateSelf(pge);
                }
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                if(!bIsChild) pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),pBackgroundColor);
                for(auto &i : vElements) {
                    i->DrawSelf(pge);
                }
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };

        public:
            void AddElement(Renderable* r, olc::vf2d Location) {
                r->vLocation = olc::vf2d(vLocation.x+Location.x,vLocation.y+1);
                r->bIsChild = true;
                vElements.push_back(r);
            }

            template <typename T>
            std::vector<T> GetAllOfClass() {
                std::vector<T> v;
                for(auto &i : vElements) {
                    if(dynamic_cast<T>(i) != nullptr) v.push_back(dynamic_cast<T>(i));
                }
                return v;
            }
        
        public:
            std::vector<Renderable*> vElements;
    };



    // Visualise a std::Map
    class TextMap : public Renderable  
    {
        public:
            // int iWidth = 100;
            // int iHeight = 16;
            int iType = 13;

            TextMap(olc::vf2d vPos) {
                vLocation = vPos;
                iHeight = 28;
            }

            TextMap() {
                vLocation = olc::vf2d(0,0);
            }

            int GetIndexFromKey(std::string s) {
                int c = 0;
                for(auto &i : vMap) {
                    if(i.first == s) return c;
                    c++;
                }
                return -1;
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                for(auto &i : vElements) {
                    i.second->UpdateSelf(pge);
                    vMap[GetIndexFromKey(i.first)].second = dynamic_cast<XOLCTS::TextInput*>(dynamic_cast<XOLCTS::CustomLine*>(i.second)->vElements[1])->sInputText;                    
                    vMap[GetIndexFromKey(i.first)].first = dynamic_cast<XOLCTS::TextInput*>(dynamic_cast<XOLCTS::CustomLine*>(i.second)->vElements[0])->sInputText;
                    i.first = dynamic_cast<XOLCTS::TextInput*>(dynamic_cast<XOLCTS::CustomLine*>(i.second)->vElements[0])->sInputText;
                }
                btNewValue->vLocation = vLocation+olc::vf2d(4,0);
                btNewValue->vLocation.y = iHeight-4;
                btNewValue->iWidth = iWidth-8;
                btNewValue->UpdateSelf(pge);
                if(btNewValue->bPressed) AddValue("Name","Color");
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                if(!bIsChild) pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),pBackgroundColor);
                for(auto &i : vElements) {
                    i.second->DrawSelf(pge);
                }
                btNewValue->DrawSelf(pge);
                
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        
        public:
            std::string sText = "Text";
            olc::Pixel pColor = olc::Pixel(128,128,128);
            // std::map<std::string,Renderable*> vElements;
            // std::map<std::string,std::string> vMap;
            std::vector<std::pair<std::string,Renderable*>> vElements;
            std::vector<std::pair<std::string,std::string>> vMap;
            Button* btNewValue = new Button(olc::vf2d(vLocation.x,0),false,"Add Value");


            void AddElement(std::string key, XOLCTS::Renderable* e) {
                e->vLocation = olc::vf2d(vLocation.x,vLocation.y+(iHeight-28)+1);
                e->bIsChild = true;
                e->iWidth = iWidth-16;
                iHeight += e->iHeight+1;
                vElements.push_back(std::pair<std::string,Renderable*>{key,e});
            };

            void AddValue(std::string key, std::string value) {
                vMap.push_back(std::pair<std::string,std::string>{key,value});
                CustomLine* c = new CustomLine();
                c->iHeight = 16;
                AddElement(key,c);
                // std::cout << "first shit" << std::endl;
                c->AddElement(new TextInput(olc::vf2d(0,0),key),olc::vf2d(2,0));
                // std::cout << "second shit" << std::endl;
                c->vElements[0]->iWidth = iWidth/2-6;
                // std::cout << "third shit" << std::endl;

                c->AddElement(new TextInput(olc::vf2d(0,0),value),olc::vf2d(iWidth/2+2,0));
                // std::cout << "fourth shit" << std::endl;
                c->vElements[1]->iWidth = iWidth/2-6;
                // std::cout << "fith shit" << std::endl;
                bRequireUpdate = true;
                iHeightChange += 16;
                // iHeight += 8;

                return;
            }
    };



    // Simple Color Picker
    class RoundColorPicker : public Renderable  
    {
        public:
            // int iWidth = 128;
            // int iHeight = 128;
            int iType = 14;

            RoundColorPicker(olc::vf2d vPos) {
                vLocation = vPos;
                iWidth = 96;
                iHeight = 96;
            }

            RoundColorPicker(olc::vf2d vPos, olc::vf2d v) {
                vLocation = vPos;
                iWidth = 96;
                iHeight = 96;
                iRadius = (int)std::min(iWidth,iHeight)/2;
                vDisplacedLocation = vLocation + olc::vf2d(iWidth/2,iHeight/2);
                vLocationOnCircle = v*iRadius;
                vTemp = v;
                pChosenColor = HSLtoRGB(180+RadToDeg(AngleFromPoints(olc::vf2d(0,0),vLocationOnCircle)),1,Map(mag(vLocationOnCircle,olc::vf2d(0,0)),0,iRadius,1,0));
                // UpdateSelf(pge);
            }

            RoundColorPicker() {
                vLocation = olc::vf2d(0,0);
                iWidth = 96;
                iHeight = 96;
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                iRadius = (int)std::min(iWidth,iHeight)/2;
                vDisplacedLocation = vLocation + olc::vf2d(iWidth/2,iHeight/2);
                if((pge->GetMouse(0).bPressed && mag(vDisplacedLocation,pge->GetMousePos()) < iRadius) || bHoldingPoint) {
                    bHoldingPoint = true;
                    vLocationOnCircle = pge->GetMousePos()-vDisplacedLocation;
                    if(mag(vDisplacedLocation,pge->GetMousePos()) > iRadius) {
                        vLocationOnCircle = EndFromAngle(olc::vf2d(0,0),180+RadToDeg(AngleFromPoints(vDisplacedLocation,pge->GetMousePos())),iRadius);
                    }
                    
                } else {
                    bHoldingPoint = false;
                }
                if(pge->GetMouse(0).bReleased) bHoldingPoint = false;
                pChosenColor = HSLtoRGB(180+RadToDeg(AngleFromPoints(olc::vf2d(0,0),vLocationOnCircle)),1,Map(mag(vLocationOnCircle,olc::vf2d(0,0)),0,iRadius,1,0));
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                if(!bIsChild) pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),pBackgroundColor);
                for (int r = 0; r < 360; r++) {
                    for (float l = 0; l < iRadius; l += 0.25) {
                        pge->Draw(EndFromAngle(vDisplacedLocation,r,l),HSLtoRGB(r,1,1.0f-(1.0f/iRadius*l)));
                    }
                }
                pge->DrawCircle(vLocationOnCircle+vDisplacedLocation,(int)iRadius/16+1,olc::WHITE);
                pge->DrawCircle(vLocationOnCircle+vDisplacedLocation,(int)iRadius/16,olc::BLACK);
                pge->DrawLine(vTemp+vDisplacedLocation,vDisplacedLocation);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        
        public:
            olc::vf2d vLocationOnCircle;
            int iRadius = 0;
            olc::vf2d vDisplacedLocation;
            olc::vf2d vTemp;
            bool bHoldingPoint = false;
            olc::Pixel pChosenColor = olc::Pixel(255,255,255);

    };



    // Basicly a renderable for shortcuts
    class Action : public Renderable  
    {
        public:
            // int iWidth = 100;
            int iHeight = 16;
            int iType = 8;

            Action(olc::vf2d vPos,KeyboardShortcut *Shortcut, std::string ActionText, int Width) {
                vLocation = vPos;
                ksShortcut = Shortcut;
                sActionName = ActionText;
                vdActionFunction = Shortcut->vdActionFunction;
                bHasShortcut = true;
                iWidth = Width;
                std::string s = "";
                for(int i = 0; i<ksShortcut->vKeys.size()-1;i++) {
                    s += GetKeyName(ksShortcut->vKeys[i])+"+";
                }
                s += GetKeyName(ksShortcut->vKeys[ksShortcut->vKeys.size()-1]);
                sShortcutText = s;
            }

            
            Action(olc::vf2d vPos,void (*Function)(), std::string ActionText, int Width) {
                vLocation = vPos;
                sActionName = ActionText;
                vdActionFunction = Function;
                iWidth = Width;
                bHasShortcut = false;
            }

            Action() {
                vLocation = olc::vf2d(0,0);
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                if(pge->GetMouse(0).bPressed && InBounds(vLocation,vLocation+olc::vf2d(iWidth,iHeight),pge->GetMousePos())) {vdActionFunction(); bUsesRenderable = true;}
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                pge->DrawString(vLocation+olc::vf2d(0,(iHeight-8)/2),sActionName,pColor);
                if(!bHasShortcut) return; 
                pge->DrawString(vLocation+olc::vf2d(iWidth-(8*sShortcutText.size())-8,(iHeight-8)/2),sShortcutText,pColor);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        
        public:
            std::string sActionName = "Action";
            std::string sShortcutText = "";
            KeyboardShortcut *ksShortcut;
            bool bHasShortcut;
            olc::Pixel pColor = olc::Pixel(128,128,128);
            void (*vdActionFunction)();

    };



    // /!\ HAVE TO REWORK /!\ Allows to create an object 
    class Object : public Renderable
    {
        public:
            bool bActive;
            int iWidth = 100;
            int iHeight = 16;
            int iType = 6;

            Object(olc::vf2d VPos, std::string name, void (*Function)(Object*), olc::PixelGameEngine* pge) {
                vLocation = VPos;
                sObjectText = name;
                vdFunction = Function;
                ptrPGE = pge;
            };

            virtual void Initialize() {
                std::cout << "Renderable Initialised." << std::endl;
                return;
            };

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                olc::HWButton MouseButttonStateArray[3] {pge->GetMouse(0),pge->GetMouse(1),pge->GetMouse(2)};
                olc::vf2d MousePos = pge->GetMousePos();
                if(bPressed) {bPressed = false;}
                if(bReleased) {bReleased = false;}
                if(!bKeepValue) {bObjectValue = false;}
                if(XOLCTS::InBounds(vLocation+olc::vf2d((iWidth-16)*!(!bKeepValue && !bShowRect),0),olc::vf2d(vLocation.x+iWidth,vLocation.y+16),MousePos) || bHeld) {
                    bHoovered = true;
                    if(MouseButttonStateArray[0].bPressed) {
                        bPressed = true;
                        bUsesRenderable = true;
                    }
                    if(MouseButttonStateArray[0].bPressed || (bHeld && !bKeepValue)) {
                        bObjectValue = !bObjectValue;
                        bHeld = true;
                    }
                } else {
                    bHoovered = false;
                }
                if(!MouseButttonStateArray[0].bHeld && bHeld) {bReleased = true;bHeld = false;}
                if(bPressed) vdFunction(this);
                if(bHeld) bUsesRenderable = true;
                if(bReleased) bUsesRenderable = false;
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                pge->SetPixelMode(olc::Pixel::ALPHA);
                pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),olc::Pixel(0,0,0,200));
                pge->DrawRect(vLocation,olc::vf2d(iWidth,iHeight),olc::Pixel(100,100,100,120));
                if(bKeepValue || bShowRect) {
                    pge->FillRect(vLocation+olc::vf2d(iWidth-16,0),olc::vf2d(15,15),olc::Pixel(0+(180*bObjectValue)+(50*bHoovered),0+(164*bObjectValue)+(50*bHoovered),0+(50*bHoovered)));
                    pge->DrawRect(vLocation+olc::vf2d(iWidth-16,0),olc::vf2d(15,15),olc::Pixel(200,200,200,150));
                    pge->DrawString(vLocation+olc::vf2d(4,4),sObjectText,olc::Pixel(255,255,255));
                }
                if(!bKeepValue && !bShowRect) {
                    pge->FillRect(vLocation,olc::vf2d(iWidth,iHeight),olc::Pixel(0+(130*bObjectValue)+(50*bHoovered),0+(114*bObjectValue)+(50*bHoovered),0+(50*bHoovered)));
                    pge->DrawRect(vLocation,olc::vf2d(iWidth-1,iHeight-1),olc::Pixel(200,200,200,150));
                    pge->DrawString(vLocation+olc::vf2d(4,4),sObjectText,olc::Pixel(50+(155*!bObjectValue),50+(155*!bObjectValue),50+(155*!bObjectValue)));
                }
                pge->SetPixelMode(olc::Pixel::NORMAL);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        public:
            bool bObjectValue = false;
            bool bKeepValue = false;
            bool bHeld = false;
            bool bPressed = false;
            bool bReleased = false;
            bool bHoovered = false;
            bool bShowRect = false;
            olc::PixelGameEngine* ptrPGE;
            void (*vdFunction)(Object*);
            // double (*dFunction)(double);
            // olc::vf2d vObjectLocation;
            // int (*CreateFunction)(double (*Function)(double), olc::vf2d Location);
            std::string sObjectText = "Object";
    };



    // Noise wraper class.
    class Noise
    {
        public:
            FastNoiseLite nNoise;

            float fFrequency = 2.5;
            float fOctaves = 5.0;
            float fPersistence = 0.5f;
            float fLacunarity = 2.5;
            float fZoom = 1;
            float z = 0;

            XOLCTS::Slider* slFrequency;
            XOLCTS::Slider* slOctaves;
            XOLCTS::Slider* slPersistence;
            XOLCTS::Slider* slLacunarity;
            XOLCTS::Slider* slz;
            XOLCTS::Slider* slZoom;

            XOLCTS::Menu mnNoiseParamitters;

            Noise() {

                nNoise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);

                mnNoiseParamitters.vLocation = olc::vf2d(0,0);
                
                mnNoiseParamitters.sMenuName = "Noise";
                
                slFrequency = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,2.5,&fFrequency,"Freq.");
                
                slOctaves = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,5,&fOctaves,"Octv.");
                
                slPersistence = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,0.5,&fPersistence,"Prst.");
                
                slLacunarity = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,2.5,&fLacunarity,"Lacu.");
                
                slz = new XOLCTS::Slider(olc::vf2d(0,0),0,360,100,0,&z,"Z");
                
                slZoom = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,1,&fZoom,"Zoom");

                mnNoiseParamitters.AddElement(slFrequency);

                mnNoiseParamitters.AddElement(slOctaves);

                mnNoiseParamitters.AddElement(slPersistence);

                mnNoiseParamitters.AddElement(slLacunarity);

                mnNoiseParamitters.AddElement(slz);		
                
                mnNoiseParamitters.AddElement(slZoom);

            }

            Noise(FastNoiseLite::NoiseType Type,bool bOverrideZSlider) {

                nNoise.SetNoiseType(Type);

                mnNoiseParamitters.vLocation = olc::vf2d(0,0);
                
                mnNoiseParamitters.sMenuName = "Noise";
                
                slFrequency = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,2.5,&fFrequency,"Freq.");
                
                slOctaves = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,5,&fOctaves,"Octv.");
                
                slPersistence = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,0.5f,&fPersistence,"Prst.");
                
                slLacunarity = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,2.5,&fLacunarity,"Lacu.");
                
                if(!bOverrideZSlider) slz = new XOLCTS::Slider(olc::vf2d(0,0),0,360,100,0,&z,"Z");
                
                slZoom = new XOLCTS::Slider(olc::vf2d(0,0),0,10,100,1,&fZoom,"Zoom");

                mnNoiseParamitters.AddElement(slFrequency);

                mnNoiseParamitters.AddElement(slOctaves);

                mnNoiseParamitters.AddElement(slPersistence);

                mnNoiseParamitters.AddElement(slLacunarity);
                
                if(!bOverrideZSlider) mnNoiseParamitters.AddElement(slz);		
                
                mnNoiseParamitters.AddElement(slZoom);

            }
            
            float GetNoise(float x,float y) {
                return XOLCTS::GetNoiseAdvanced(nNoise,x/fZoom,y/fZoom,z,fFrequency,fOctaves,fPersistence,fLacunarity);
            }

            float GetNoise(float x,float y, float oz) {
                return XOLCTS::GetNoiseAdvanced(nNoise,x/fZoom,y/fZoom,oz,fFrequency,fOctaves,fPersistence,fLacunarity);
            }
  

    };



    // Particle class thing
    class Particle {
        public:
            olc::vf2d vLocation;
            olc::vf2d vVelocity;
            olc::vf2d vAcceleration;
            olc::vf2d vOffset;

            bool bDead = false;
            bool bVelDecreases = true;
            bool bRotates;
            bool bChangesHue;
            bool bDisapears = false;
            int iOriginalLifetime = -1;
            int iLifetime = -1;


            float fAirResistance = 0.95;
            float fRotationDegree = 0;
            float fHue = 0;
            float fHueChangeSpeed = 1;

            olc::Decal* decParticle;


            Particle(olc::vf2d loc, olc::vf2d vel, olc::vf2d acc, int l, olc::Decal* d, bool rotates, bool huechange) {
                vLocation = loc;
                vVelocity = vel;
                vAcceleration = acc;
                iOriginalLifetime = l;
                iLifetime = l;
                decParticle = d;
                bRotates = rotates;
                bChangesHue = huechange;
                fHueChangeSpeed = 1;
                return;
            }

            void Create(olc::vf2d loc, olc::vf2d vel, olc::vf2d acc, int l, olc::Decal* d, bool rotates, bool huechange) {
                vLocation = loc;
                vVelocity = vel;
                vAcceleration = acc;
                iOriginalLifetime = l;
                iLifetime = l;
                decParticle = d;
                bRotates = rotates;
                bChangesHue = huechange;
                bDead = false;
                fRotationDegree = 0;
                fHue = 0;
                fHueChangeSpeed = 1;
                return;
            }

            void SetLifetime(int i) {
                iOriginalLifetime = i;
                iLifetime = i;
                return;
            }

            void AddAcceleration(olc::vf2d Acc) {
                vAcceleration += Acc;
                return;
            }

            void Update(olc::PixelGameEngine *pge) {
                if(bDead) return;
                vVelocity += vAcceleration;
                vAcceleration = {0,0};
                vLocation += vVelocity;
                vVelocity *= fAirResistance*(1*bVelDecreases);
                if(bRotates) fRotationDegree = (fRotationDegree > 360) ? 0 : fRotationDegree+0.125;
                if(bChangesHue && fHueChangeSpeed > 0) fHue = (fHue > 360) ? 0 : fHue + fHueChangeSpeed;
                if(bChangesHue && fHueChangeSpeed < 0) fHue = (fHue < 0) ? 360 : fHue + fHueChangeSpeed;
                if(!InBounds({0,0},{(float)pge->ScreenWidth(),(float)pge->ScreenHeight()},vLocation)) bDead = true;
                if(iLifetime == 0) bDead = true;
                else if(iLifetime != -1) iLifetime--;
                return;
            }

            void Show(olc::PixelGameEngine *pge) {
                if(bDead) return;
                olc::Pixel pCol = (bChangesHue) ? HSLtoRGB(fHue,0.5,0.5) : olc::WHITE;
                pCol.a = Map(iLifetime,0,iOriginalLifetime,0,255);
                pge->DrawRotatedDecal(vLocation+vOffset+vShake,decParticle,fRotationDegree,{(float)decParticle->sprite->width/2,(float)decParticle->sprite->height/2}, vZoom, pCol);
                return;
            }

            Particle* CreateParticle(olc::vf2d loc, olc::vf2d vel, olc::vf2d acc, int l, olc::Decal* d, bool rotates, bool huechange);
    };



    // Simple Console Inplementation
    class Console : public Renderable  
    {
        public:
            // int iWidth = 100;
            // int iHeight = 16;
            int iType = 15;

            Console(olc::vf2d vPos,olc::vf2d vSize,olc::PixelGameEngine* pge) {
                ptrPGE = pge;
                vLocation = vPos;
                bEditable = true;
                Resize(vSize);
                GenerateSprite(pge);
                RegisterBaseCommands();
            } 

            Console(olc::vf2d vPos,float x, float y, olc::PixelGameEngine* pge) {
                ptrPGE = pge;
                vLocation = vPos;
                bEditable = true;
                Resize(x,y);
                GenerateSprite(pge);
                RegisterBaseCommands();
            } 


            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                bUpperCase = GetCapsLockState();
                if(pge->GetKey(olc::SHIFT).bHeld) bUpperCase = true;
                iMaxCharsIndex = floor(iWidth/8);
                if(!bEditable) return;
                if(pge->GetMouse(0).bPressed) {
                    bSelected = InBounds(vLocation,vLocation+olc::vf2d(iWidth,iHeight),pge->GetMousePos());
                    iCursorIndex = (int)((pge->GetMousePos().x-vLocation.x)/8)*bSelected;
                }
                while(vConsoleBuffer.size() > floor(iHeight/16)-1) vConsoleBuffer.erase(vConsoleBuffer.begin());
                if(bSelected) {
                    bUsesRenderable = true;

                    iCursorIndex = iCursorIndex+(1*pge->GetKey(olc::RIGHT).bPressed || (pge->GetKey(olc::RIGHT).bHeld && iKeyTimeout == 0));
                    iCursorIndex = iCursorIndex-(1*pge->GetKey(olc::LEFT).bPressed || (pge->GetKey(olc::LEFT).bHeld && iKeyTimeout == 0));
                    if(pge->GetKey(olc::LEFT).bHeld && iKeyTimeout == 0) iKeyTimeout = 10;
                    if(pge->GetKey(olc::RIGHT).bHeld && iKeyTimeout == 0) iKeyTimeout = 10;



                    if(iCursorIndex < 0) iCursorIndex = 0;
                    if(iCursorIndex > sInputText.size()) iCursorIndex = sInputText.size();
                    if(iCursorIndex-iDisplayOffset > iMaxCharsIndex) iDisplayOffset++;
                    if(iCursorIndex-iDisplayOffset < 1) iDisplayOffset = (iDisplayOffset > iMaxCharsIndex-1) ? iDisplayOffset-iMaxCharsIndex : 0;



                    std::unordered_set<olc::Key> vKeysPressed;
                    for(std::map<size_t, uint8_t>::iterator it = olc::mapKeys.begin(); it != olc::mapKeys.end(); ++it) {
                        if(pge->GetKey(olc::Key(it->second)).bPressed) {vKeysPressed.insert(olc::Key(it->second));}
                    }



                    std::vector<std::string> sUserInput = GetKeys(vKeysPressed,bUpperCase);


                    if(pge->GetKey(olc::UP).bPressed && iCommandHistoryIndex > 0) {iCommandHistoryIndex--; LoadCommandFromHistory(&sUserInput);}
                    if(pge->GetKey(olc::DOWN).bPressed && iCommandHistoryIndex < vLastCommandsIssued.size()) {iCommandHistoryIndex++; LoadCommandFromHistory(&sUserInput);}


                    if(pge->GetKey(olc::CTRL).bHeld && pge->GetKey(olc::V).bPressed) {sUserInput.clear(); sUserInput = XOLCTS::SplitString(GetClipboard(),"\n");}
                    if(pge->GetKey(olc::SHIFT).bHeld && pge->GetKey(olc::Key(35)).bPressed) {sUserInput.clear(); sUserInput.push_back("(");}
                    if(pge->GetKey(olc::SHIFT).bHeld && pge->GetKey(olc::Key(36)).bPressed) {sUserInput.clear(); sUserInput.push_back(")");}
                    if(pge->GetKey(olc::SHIFT).bHeld && pge->GetKey(olc::Key(34)).bPressed) {sUserInput.clear(); sUserInput.push_back("\\");}


                    for(auto &i : sUserInput) {
                        sInputText.insert(iCursorIndex,i);
                        iCursorIndex += i.size(); 
                        iDisplayOffset += i.size();
                    }
                    iKeyTimeout -= (iKeyTimeout < 1) ? 0 : 1*(pge->GetElapsedTime()*2);
                    if(((pge->GetKey(olc::BACK).bHeld && iKeyTimeout == 0) || pge->GetKey(olc::BACK).bPressed) && iCursorIndex > 0 ) {sInputText.erase(sInputText.begin()+iCursorIndex-1); iCursorIndex--; iKeyTimeout = iKeyTimeoutDuration*2*(pge->GetKey(olc::BACK).bPressed+1);} 
                    if(((pge->GetKey(olc::DEL).bHeld && iKeyTimeout == 0) || pge->GetKey(olc::DEL).bPressed) && iCursorIndex < sInputText.size() && iKeyTimeout == 0) {sInputText.erase(sInputText.begin()+iCursorIndex); iKeyTimeout = iKeyTimeoutDuration*2*(pge->GetKey(olc::DEL).bPressed+1);} 
                    bReturned = pge->GetKey(olc::ENTER).bPressed; 
                    if(bReturned) CheckCommands();
                }
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                std::string sDrawnString = sInputText;
                if(sInputText.size() > iMaxCharsIndex) {
                    sDrawnString = "";
                    for(int i = iDisplayOffset*bSelected; (i < (iDisplayOffset*bSelected)+iMaxCharsIndex && i < sInputText.size()); i++) {
                        sDrawnString += sInputText[i];
                    }
                } else iDisplayOffset = 0;


                pge->DrawDecal(vLocation,sConsoleDecal);
                int igd = 0;
                if(!bDebug) {
                    for(int i = vConsoleBuffer.size()-1; i > -1; i--) {
                        if(vConsoleBuffer[(vConsoleBuffer.size()-1)-i].substr(0,5) == "~<d>~") {
                            igd++;
                        }
                    }
                }

                int ih = vConsoleBuffer.size()-1;
                for(int i = vConsoleBuffer.size()-1; i > -1; i--) {
                    if(vConsoleBuffer[(vConsoleBuffer.size()-1)-i].substr(0,5) == "~<d>~") {
                        if(bDebug){
                             pge->DrawStringDecal(vLocation+olc::vf2d(4,iHeight-(16*(ih-igd+1))-12),vConsoleBuffer[(vConsoleBuffer.size()-1)-i].substr(5,vConsoleBuffer[(vConsoleBuffer.size()-1)-i].size()-5),pTextColor);
                        } else ih++;
                        ih--;
                    } else {
                        pge->DrawStringDecal(vLocation+olc::vf2d(4,iHeight-(16*(ih-igd+1))-12),vConsoleBuffer[(vConsoleBuffer.size()-1)-i],pTextColor);
                        ih--;
                    }
                }
            

                if(sInputText.size() > 0) {
                    pge->DrawStringDecal(vLocation+olc::vf2d(0,iHeight-16)+olc::vf2d(4,4),sDrawnString,pTextColor);
                } else {
                    sHintText = std::string(">> Type Command").substr(0,iMaxCharsIndex);
                    pge->DrawStringDecal(vLocation+olc::vf2d(0,iHeight-16)+olc::vf2d(4,4),sHintText,pTextColor*0.75);
                }
                if(bSelected) pge->FillRectDecal(vLocation+olc::vf2d(0,iHeight-16)+olc::vf2d(4+(8*(iCursorIndex-iDisplayOffset)),2),olc::vf2d(1,iHeight-2),olc::WHITE);
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };

        private:
            bool bUpperCase = false;
        
        public:


            std::string sHintText = "";
            std::string sInputText = "";
            bool bSelected = false;
            int iCursorIndex = 0;
            int iDisplayOffset = 0;
            int iMaxCharsIndex = 0;
            int iKeyTimeout = 0;
            int iKeyTimeoutDuration = 3;
            int iCommandHistoryIndex = 0;
            bool bEditable;
            bool bDebug = false;
            olc::Pixel pTextColor = olc::Pixel(128,128,128);
            bool bReturned = false;
            std::map<std::string, void (*)(XOLCTS::Console*)> vCommands;
            std::map<std::string, void (XOLCTS::Console::*)(XOLCTS::Console*)> vInternalCommands;
            std::vector<std::string> vLastCommandsIssued;
            std::vector<std::string> vConsoleBuffer;
            std::vector<std::string> vArgs;
            olc::PixelGameEngine* ptrPGE;
            olc::Decal* sConsoleDecal;
            olc::Sprite* sConsoleSprite;

            void RegisterBaseCommands() {
                RegisterCommand("ToggleDebug",conToggleDebug);
                RegisterCommand("Clear",conClearConsole);
                return;
            }

            static void conToggleDebug(XOLCTS::Console* _con) {
                _con->bDebug = !_con->bDebug;
                if(_con->bDebug) _con->WriteText("Enabled debug output.");
                else _con->WriteText("Disabled debug output.");
            }

            static void conClearConsole(XOLCTS::Console* _con) {
                _con->vConsoleBuffer.clear();
            }

            void Resize(float x, float y) {
                iWidth = x;
                iHeight = y;
                sConsoleSprite = new olc::Sprite(iWidth,iHeight);
                GenerateSprite(ptrPGE);
                sConsoleDecal = new olc::Decal(sConsoleSprite);
            };
            
            void Resize(olc::vf2d size) {
                iWidth = size.x;
                iHeight = size.y;
                sConsoleSprite = new olc::Sprite(iWidth,iHeight);
                GenerateSprite(ptrPGE);
                sConsoleDecal = new olc::Decal(sConsoleSprite);
            };

            void GenerateSprite(olc::PixelGameEngine* pge) {
                pge->SetDrawTarget(sConsoleSprite);
                pge->FillRect(olc::vf2d(0,0),olc::vf2d(iWidth,iHeight),olc::Pixel(pBackgroundColor.r*(!bEditable+1),pBackgroundColor.g*(!bEditable+1),pBackgroundColor.b*(!bEditable+1),191));
                pge->DrawRect(olc::vf2d(0,0),olc::vf2d(iWidth,iHeight),olc::BLACK);
                pge->DrawRect(olc::vf2d(1,1),olc::vf2d(iWidth-1,iHeight-1),pBackgroundColor*0.75);
                pge->FillRect(olc::vf2d(0,iHeight-16),olc::vf2d(iWidth,iHeight),olc::Pixel(pBackgroundColor.r*(!bEditable+1),pBackgroundColor.g*(!bEditable+1),pBackgroundColor.b*(!bEditable+1),223));
                pge->DrawRect(olc::vf2d(0,iHeight-16),olc::vf2d(iWidth,iHeight),olc::Pixel(0,0,0));
                pge->DrawRect(olc::vf2d(1,iHeight-15),olc::vf2d(iWidth-1,iHeight-1),pBackgroundColor*0.75);
                pge->SetDrawTarget(nullptr);
                return;
            }

            void LoadCommandFromHistory(std::vector<std::string>* sUserInput) {
                sUserInput->clear(); 
                if(iCommandHistoryIndex < vLastCommandsIssued.size() && iCommandHistoryIndex > -1) {
                    sUserInput->push_back(vLastCommandsIssued[iCommandHistoryIndex]);
                } else sUserInput->push_back("");
                sInputText = "";
                iDisplayOffset = 0;
                iCursorIndex = 0;
                return;
            }

            void WriteText(std::string s) {
                time_t now = time(0);
                tm *ltm = localtime(&now);
                std::string time = 
                "["+std::to_string(ltm->tm_hour)+
				":"+std::to_string(ltm->tm_min)+
				":"+std::to_string(ltm->tm_sec)+"] ";
                s = time+s;
                for(int i = 0; i < Constrain(floor(s.size()/((iWidth-8)/8-1))+1,1,INFINITY); i++) {
                    vConsoleBuffer.push_back(s.substr(i*((iWidth-8)/8-1),((iWidth-8)/8)-1));
                }
                vConsoleBuffer[vConsoleBuffer.size()-1] += std::string("\n");
                return;
            }

            void WriteDebug(std::string s) {
                time_t now = time(0);
                tm *ltm = localtime(&now);
                std::string time = 
                "["+std::to_string(ltm->tm_hour)+
				":"+std::to_string(ltm->tm_min)+
				":"+std::to_string(ltm->tm_sec)+"] ";
                s = time+s;
                for(int i = 0; i < Constrain(floor(s.size()/((iWidth-8)/8-1))+1,1,INFINITY); i++) {
                    vConsoleBuffer.push_back("~<d>~"+s.substr(i*((iWidth-8)/8-1),((iWidth-8)/8)-1));
                }
                vConsoleBuffer[vConsoleBuffer.size()-1] += std::string("\n");
                return;
            }

            void RegisterCommand(std::string sCommandAlias, void (*vFunction)(XOLCTS::Console*)) {
                vCommands[sCommandAlias] = vFunction;
            } 
            // void RegisterCommand(std::string sCommandAlias, void (XOLCTS::Console::*vFunction)(XOLCTS::Console*)) {
            //     vInternalCommands[sCommandAlias] = vFunction;
            // }

            void CheckCommands() {
                WriteText(">> "+sInputText);
                vArgs.clear();
                vLastCommandsIssued.push_back(sInputText);
                iCommandHistoryIndex = vLastCommandsIssued.size();
                std::string sCommand = "";
                if(sInputText.find('(') != std::string::npos && sInputText.find(')') != std::string::npos) {
                    vArgs.push_back(sInputText.substr(sInputText.find_first_of("(")+1,sInputText.size()-sInputText.find_first_of("(")-2));
                    sCommand = XOLCTS::SplitString(sInputText,"(")[0];
                } else {
                    sCommand = sInputText;
                }
                if(vCommands.count(sCommand) > 0) {
                    vCommands[sCommand](this);
                }
                vArgs.clear();
                sInputText = "";
            }
    };



    // Basic menu class. Used to list other Renderables.
    class ClickMenu : public Renderable 
    {        
        public:
            // bool bActive;
            int iWidth = 0;
            int iHeight = 16;
            int iType = 16;


            ClickMenu() {
                return;
            };

            ClickMenu(olc::vf2d Location, int Width) {
                iWidth = Width;
                vLocation = Location;
            };

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                if(!bActive) return;
                if(bOpened) { 
                    bUsesRenderable = true;
                }
                if(pge->GetMouse(1).bPressed) {
                    bOpened = true;
                    vLocation = pge->GetMousePos();
                    iMenuHeight = 2;
                    for(auto &e : vElements) {
                        e->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iMenuHeight+1);
                        e->bIsChild = true;
                        e->iWidth = iWidth-16;
                        iMenuHeight += e->iHeight+1;
                    }
                }



                if(bOpened) {
                    int index = 0;
                    iHoveredElement = -1;
                    for(auto &i : vElements) {
                        if(!i->bUpdateWithParent || !i->bActive) continue;
                        i->UpdateSelf(pge);
                        if(i->bRequireUpdate) {iMenuHeight += i->iHeightChange; i->bRequireUpdate = false; i->iHeightChange = 0;};
                        if(i->bSelectable && InBounds(i->vLocation-olc::vf2d(8,0),i->vLocation+olc::vf2d(iWidth,i->iHeight),pge->GetMousePos())) iHoveredElement = index;
                        index++;
                    }
                }
                if(pge->GetMouse(0).bReleased && (!InBounds(vLocation,vLocation+olc::vf2d(iWidth,iHeight),pge->GetMousePos()) || iHoveredElement != -1 && dynamic_cast<Action*>(vElements[iHoveredElement]) != nullptr)) {
                    bOpened = false;
                }

                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                if(!bActive || !bOpened) return; // Stops drawing if menu isn't opened (prevents from drawing the menu elements).
                pge->SetPixelMode(olc::Pixel::ALPHA);
                pge->FillRect(vLocation,olc::vf2d(iWidth,iMenuHeight+3),pBackgroundColor); // Drawing the menu text box.
                pge->DrawRect(vLocation,olc::vf2d(iWidth,iMenuHeight+3),olc::BLACK); // Draws menu box outline.
                if(iHoveredElement != -1) pge->FillRect(vElements[iHoveredElement]->vLocation-olc::vf2d(7,0),olc::vf2d(iWidth-2,vElements[iHoveredElement]->iHeight),olc::Pixel(50, 50, 50)); // Draws element highlight if necessary.

                pge->SetPixelMode(olc::Pixel::NORMAL);
                for(auto &i : vElements) {
                    i->DrawSelf(pge);
                    pge->SetPixelMode(olc::Pixel::ALPHA);
                    if(!i->bActive) pge->FillRect(i->vLocation,olc::vf2d(i->iWidth,i->iHeight),olc::Pixel(pBackgroundColor.r/2,pBackgroundColor.g/2,pBackgroundColor.b/2,223));
                    pge->SetPixelMode(olc::Pixel::NORMAL);
                }
                return;
            };

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };


        public:
            std::vector<XOLCTS::Renderable*> vElements;
            olc::Pixel pColor = olc::Pixel(37,37,37);
            int iMenuHeight = 2;
            int iSubMenuTimer = 0;
            int iHoveredElement = -1; 
            bool bOpened = false;

            void AddElement(XOLCTS::Renderable* e) {
                e->vLocation = olc::vf2d(vLocation.x+8,vLocation.y+iHeight+iMenuHeight+1);
                e->bIsChild = true;
                e->iWidth = iWidth-16;
                iMenuHeight += e->iHeight+1;
                vElements.push_back(e);
            };

            void AddSeparator() {
                AddElement(new Separator(olc::vf2d(0,0),iWidth-16));
            }

            SubMenu* AddSubmenu(olc::vf2d Location, int Width, int MenuWidth, std::string Text, int MenuID) {
                SubMenu* i = new SubMenu(Location, Width, MenuWidth, Text, MenuID);
                AddElement(i);
                return  i;
            }

            template <typename T>
            std::vector<T> GetAllOfClass() {
                std::vector<T> v;
                for(auto &i : vElements) {
                    if(dynamic_cast<T>(i) != nullptr) v.push_back(dynamic_cast<T>(i));
                }
                return v;
            }

            template <typename T>
            std::vector<int> GetAllIndexOfClass() {
                int c = 0;
                std::vector<int> v;
                for(auto &i : vElements) {
                    if(dynamic_cast<T>(i) != nullptr) v.push_back(c);
                    c++;
                }
                return v;
            }

            int GetIndexOfItem(Renderable* r) {
                int c = 0;
                for(auto &i : vElements) {
                    if(r == i) return c;
                    c++;
                }
                return -1;
            }

            void ResetMenu() {
                for(auto &i : vElements) {
                    delete i;
                }
                vElements.clear();
                
                iMenuHeight = 2;
            }
    };



    // Controller Base class for controlable objects
    class Controller {
        public:
            virtual void Update(float fElapsedTime) {
                return;
            }


        public:
            olc::PixelGameEngine* ptrPGE;
            double dZoom;
            olc::vf2d vLocation;		 

    };



    // Default mous controller for map like objects
    class DefaultMouseController : public Controller {
        public:
            DefaultMouseController(olc::PixelGameEngine* pge) {
                ptrPGE = pge;
                vLocation = olc::vf2d(0,0);
                dZoom = 1;
            }

        
        public:
            virtual void Update(float fElapsedTime) {
                if(ptrPGE->GetMouse(0).bPressed) vClickMouseLoc = olc::vf2d(ptrPGE->GetMousePos());
                if(ptrPGE->GetMouse(0).bHeld) {vLocation -= (vClickMouseLoc-ptrPGE->GetMousePos())/dZoom; vClickMouseLoc = olc::vd2d(ptrPGE->GetMousePos());}
                olc::vf2d v1 = ScreenToWorld(ptrPGE->GetMousePos());
                if((double)ptrPGE->GetMouseWheel() > 1) dZoom /= 0.90;
                if((double)ptrPGE->GetMouseWheel() < -1) dZoom *= 0.90;
                olc::vf2d v2 = ScreenToWorld(ptrPGE->GetMousePos());
                vLocation -= v1-v2;
            }


        public:
            olc::vf2d vClickMouseLoc;
            
            olc::vf2d WorldToScreen(olc::vf2d vWorldCoord) {
                return ((vWorldCoord)-(vLocation))*dZoom;
            }

            olc::vf2d ScreenToWorld(olc::vf2d vScreenCoord) {
                return (vScreenCoord/dZoom+(vLocation));
            }

    };



    class DrawingUnit {
        public:
            DrawingUnit() {
                return;
            }

            virtual void Draw(olc::vf2d vLocation,int v, olc::PixelGameEngine* pge, double Size) {
                return;
            };


        public:
            olc::PixelGameEngine* ptrPGE;
    };



    // Simple Texture array 
    class TextureArray : public DrawingUnit {
        public:

            TextureArray() {
                return;
            }

            void Draw(olc::vf2d vLocation,int v, olc::PixelGameEngine* pge, double Size) {
                pge->DrawSprite(vLocation,vcSprites[v],Size);
                return;
            }

            olc::Sprite* NewSprite() {
                vcSprites.push_back(new olc::Sprite(16,16));
                return vcSprites[vcSprites.size()-1];
            }

        public:
            std::vector<olc::Sprite*> vcSprites;
    };



    // Simple Decal array 
    class DecalArray : public DrawingUnit {
        public:

            DecalArray() {
                return;
            }

            void Draw(olc::vf2d vLocation,int v, olc::PixelGameEngine* pge, double Size) {
                pge->DrawDecal(vLocation,vcDecal[v].get(),olc::vf2d(Size,Size));
                // pge->DrawRect(vLocation,olc::vf2d(16*Size,16*Size),olc::WHITE);
                return;
            }

            olc::Decal* NewDecal() {
                vcSprites.push_back(std::make_unique<olc::Sprite>(16,16));
                vcDecal.push_back(std::make_unique<olc::Decal>(vcSprites[vcSprites.size()-1].get()));
                return vcDecal[vcDecal.size()-1].get();
            }
            
            olc::Decal* NewDecalFromFile(std::string s) {
                vcSprites.push_back(std::make_unique<olc::Sprite>(16,16));
                vcSprites[vcSprites.size()-1].get()->LoadFromFile(s);
                vcDecal.push_back(std::make_unique<olc::Decal>(vcSprites[vcSprites.size()-1].get()));
                return vcDecal[vcDecal.size()-1].get();
            }

        public:
            std::vector<std::unique_ptr<olc::Sprite>> vcSprites;
            std::vector<std::unique_ptr<olc::Decal>> vcDecal;
    };



    // Basic grid implementation
    class Grid : public Renderable {

        public:

            Grid(olc::vf2d vPos,olc::vf2d TileSize, olc::vi2d Size, olc::PixelGameEngine* pge) {
                vLocation = vPos;
                vGridSize = Size;
                vTileSize = TileSize;
                vGridOffset = olc::vf2d(0,0);
                pgePTR = pge;
                cGridController = new DefaultMouseController(pgePTR);
                vcTiles.reserve(vGridSize.x);
                for(int x = 0; x < vGridSize.x; x++) {
                    vcTiles[x].reserve(vGridSize.y);
                }
                for(	int x = 0; 	x < vGridSize.x*vTileSize.x; x += vTileSize.x) {
                    for(int y = 0; 	y < vGridSize.y*vTileSize.y; y += vTileSize.y) {
                        vcTiles[(int)x/vTileSize.x][(int)y/vTileSize.y] = -1;
                    }
                }
            }

            virtual void UpdateSelf(olc::PixelGameEngine *pge) {
                cGridController->Update(pge->GetElapsedTime());
                vGridOffset = -cGridController->vLocation;
                dZoom = cGridController->dZoom;
                return;
            };

            virtual void DrawSelf(olc::PixelGameEngine *pge) {
                for(	int x =  	Constrain(ScreenToWorld(olc::vf2d(0,0)).x-((int)vGridOffset.x%(int)vTileSize.x),0,INFINITY); 	x < Constrain(ScreenToWorld(olc::vf2d(0,0)).x+((double)pge->GetWindowSize().x)/(vTileSize.x*dZoom)*vTileSize.x,0,vGridSize.x*vTileSize.x); x += vTileSize.x) {
                    for(int y =		Constrain(ScreenToWorld(olc::vf2d(0,0)).y-((int)vGridOffset.y%(int)vTileSize.y),0,INFINITY); 	y < Constrain(ScreenToWorld(olc::vf2d(0,0)).y+((double)pge->GetWindowSize().y)/(vTileSize.y*dZoom)*vTileSize.y,0,vGridSize.y*vTileSize.y); y += vTileSize.y) {
                        if(!XOLCTS::InBounds(WorldToScreen(olc::vf2d(0,0)),WorldToScreen(vGridSize*vTileSize),WorldToScreen(olc::vf2d(x,y)))) continue;
                        duGridDrawingUnit->Draw(WorldToScreen(olc::vf2d(x,y)),vcTiles[(int)x/vTileSize.x][(int)y/vTileSize.y],pge,dZoom);
                    }
                }
                return;
            };

            

        public:

            virtual int GetType() {
                return iType;
            };

            virtual int GetWidth() {
                return iWidth;
            };

            virtual int GetHeight() {
                return iHeight;
            };
        private:
            olc::vf2d vGridOffset;
            double dZoom = 1;
        public:
            olc::vi2d vGridSize;
            olc::vf2d vTileSize;
            std::vector<std::vector<int>> vcTiles;
            
            olc::PixelGameEngine* pgePTR;
            Controller* cGridController;
            DrawingUnit* duGridDrawingUnit;
            


            olc::vf2d WorldToScreen(olc::vf2d vWorldCoord) {
                return ((vWorldCoord)-(vLocation+vGridOffset))*dZoom;
            }

            olc::vf2d ScreenToWorld(olc::vf2d vScreenCoord) {
                return (vScreenCoord/dZoom+(vLocation+vGridOffset));
            }

            olc::vf2d ScreenToTile(olc::vf2d vScreenCoord) {
                return (olc::vi2d)((vScreenCoord/dZoom+(vLocation+vGridOffset))/vTileSize)*vTileSize;
            }

            int GetTileFromScreen(olc::vf2d vScreenCoord) {
                olc::vi2d v = ((vScreenCoord/dZoom+(vLocation+vGridOffset))/vTileSize);
                if(!XOLCTS::InBounds(olc::vf2d(0,0),vGridSize,v)) throw "Error : Not a valid tile!";
                return vcTiles[v.x][v.y];
            }

            void SetTileFromScreen(olc::vf2d vScreenCoord, int value) {
                olc::vi2d v = ((vScreenCoord/dZoom+(vLocation+vGridOffset))/vTileSize);
                if(!XOLCTS::InBounds(olc::vf2d(0,0),vGridSize,v)) return;
                vcTiles[v.x][v.y] = value; 
            }
            
            void SetTile(olc::vi2d vTile, int value) {
                if(!XOLCTS::InBounds(olc::vf2d(0,0),vGridSize,vTile)) return;
                vcTiles[vTile.x][vTile.y] = value;
            }

            int GetTile(olc::vi2d vTile) {
                if(!XOLCTS::InBounds(olc::vf2d(0,0),vGridSize,vTile)) throw "Error : Not a valid tile!";
                return vcTiles[vTile.x][vTile.y];
            }
    };



    #pragma endregion Classes

    

    #pragma region ClassDependentFunctions



    // Top Menu Array
    std::vector<TopMenu*> vTopMenuArray;



    // Null Menu
    TopMenu* NullMenu = new TopMenu(olc::vf2d(-1,-1), 0, 0, "Null", -1);



    // Creates a new Top Menu
    TopMenu* AddTopMenu(olc::vf2d Location, int Width, int MenuWidth, std::string Text, int MenuID) {
        TopMenu* i = new TopMenu(Location,Width,MenuWidth, Text, MenuID);
        vTopMenuArray.push_back(i);
        return i;
    }



    // Get Top Menu by ID
    TopMenu* GetTopMenuFromID(int ID) {
        for(auto &i : vTopMenuArray) {
            if(i->iMenuID == ID) {
                return i;
            }
        }
        return NullMenu;
    }



    // Updates the Top Menu
    void UpdateTopMenu(olc::PixelGameEngine* pge) {
        if(pge->GetKey(olc::Key::F1).bPressed) bTopMenuShown = !bTopMenuShown;
        if(!bTopMenuShown) {iMenuOpened = -1; return;}
        for(auto &i : vTopMenuArray) {
            i->UpdateSelf(pge);
        }
        return;
    }



    // Draws the Top Menu
    void DrawTopMenu(olc::PixelGameEngine* pge) {
        if(!bTopMenuShown) return;
        pge->FillRect(olc::vf2d(0,0),olc::vf2d(pge->ScreenWidth(),16),olc::Pixel(37,37,37));
        for(auto &i : vTopMenuArray) {
            i->DrawSelf(pge);
        }
        return;
    }



    // Create a Move Request
    void MoveObject(XOLCTS::Renderable* e, olc::vf2d vDestination, int Time, float fElapsedTime=1) {
        vcMoveRequestList.push_back(new MoveRequest{Time,e->vLocation,vDestination,(e->vLocation-vDestination)/Time,(int)vcMoveRequestList.size(),e});
        // std::cout << (int)vcMoveRequestList.size()-1 << std::endl;
    }



    // Updates all the Move Requests
    void UpdateMoveRequests() {
        // std::cout << "1" << std::endl;
        std::vector<int> vcFinishedTasks;
        for(auto i : vcMoveRequestList) {
            i->e->vLocation -= i->vIncrement;
            i->iTime--;
            if(i->iTime < 0) {
                vcFinishedTasks.push_back(i->iID);
                i->e->vLocation = i->vDestination;
            }
        }
        // std::cout << "2" << std::endl;
        int EreasedCount = 0;
        for(int i : vcFinishedTasks) {
            vcMoveRequestList.erase(vcMoveRequestList.begin()+(i-EreasedCount));
            // std::cout << "Move Request with ID "+std::to_string(i)+" ended." << std::endl;
            EreasedCount++;
        }
    }



    // Particle Vector
    std::vector<Particle*> vParticles;



    Particle* CreateParticle(olc::vf2d loc, olc::vf2d vel, olc::vf2d acc, int l, olc::Decal* d, bool rotates, bool huechange) {
        for(auto &i : vParticles) {
            if(i->bDead) {i->Create(loc, vel, acc, l, d, rotates, huechange); return i;}
        }
        vParticles.push_back(new Particle(loc, vel, acc, l, d, rotates, huechange));
        return vParticles[vParticles.size()-1];
    }



    Particle* CreateParticle(olc::vf2d loc, olc::vf2d vel, olc::vf2d acc, int l, olc::Decal* d, bool rotates, bool huechange, float huechangespeed) {
        for(auto &i : vParticles) {
            if(i->bDead) {i->Create(loc, vel, acc, l, d, rotates, huechange);i->fHueChangeSpeed = huechangespeed; return i;}
        }
        vParticles.push_back(new Particle(loc, vel, acc, l, d, rotates, huechange));
        vParticles[vParticles.size()-1]->fHueChangeSpeed = huechangespeed;
        return vParticles[vParticles.size()-1];
    }



	// Updates and Shows the particles
	void UpdateAndDrawParticles(olc::PixelGameEngine *pge) {
		for(auto &i : XOLCTS::vParticles) {
			if(i->bDead) continue;
			i->Update(pge);
			i->Show(pge);
		}
	}



	// Updates the particles
	void UpdateParticles(olc::PixelGameEngine *pge) {
		for(auto &i : XOLCTS::vParticles) {
			if(i->bDead) continue;
			i->Update(pge);
		}
	}



	// Shows the particles
	void ShowParticles(olc::PixelGameEngine *pge) {
		for(auto &i : XOLCTS::vParticles) {
			if(i->bDead) continue;
			i->Show(pge);
		}
	}



    #pragma endregion ClassDependentFunctions


}