#ifndef IMG2SCR_H
#define IMG2SCR_H
class ImgScr{
    float ca,sa;
    float m_sOx, m_sOy;
public:
    float m_tOx, m_tOy;
    float m_scale;
    float m_angle;
    int m_sWidth, m_sHeight;
    int m_tWidth, m_tHeight;

    ImgScr(int sWidth, int sHeight):m_tOx(0.f),m_tOy(0.f),m_scale(1.f),m_angle(0.f){
      m_sWidth  = sWidth;
      m_sHeight = sHeight;
      ca = 1.f; sa = 0.f;
      SetTexDims(1, 1);
      compute_sO();
    }
    ImgScr(const ImgScr& c){
        m_tOx = c.m_tOx; m_tOy = c.m_tOy;
        m_sOx = c.m_sOx; m_sOy = c.m_sOy;
        ca = c.ca; sa = c.sa;
        m_scale = c.m_scale;
        m_angle = c.m_angle;
        m_sWidth = c.m_sWidth; m_sHeight = c.m_sHeight;
        m_tWidth = c.m_tWidth; m_tHeight = c.m_tHeight;
    }

    void SetScreenDims(int sWidth, int sHeight){
        m_sWidth  = sWidth;
        m_sHeight = sHeight;
    }

    void SetTexDims(int tWidth, int tHeight){
        m_tWidth  = tWidth;
        m_tHeight = tHeight;
    }

    //transforms corrdinates (sx,sy) in screen space into (tx,ty) in texture image space
    void S2T(float sx, float sy, float& tx, float& ty){
        tx = 1.f/m_scale *( ca * sx + sa * sy) + m_sOx;
        ty = 1.f/m_scale *(-sa * sx + ca * sy) + m_sOy;
    }

    //transforms corrdinates (tx,ty) in texture image space into (tx,ty) in screen space
    void T2S(float tx, float ty, float& sx, float& sy){
        sx = m_scale * (ca * tx - sa * ty) + m_tOx;
        sy = m_scale * (sa * tx + ca * ty) + m_tOy;
    }

    //transforms corrdinates (tx,ty) in texture image space into (tx,ty) in screen space then scales it for OpenGL
    //That is lower left corner at (-1,-1) and upper left corner at (1,1)
    void T2S1(float tx, float ty, float& sx, float& sy){
        T2S(tx,ty,sx,sy);
        sx = 2.f / float(m_sWidth)  * sx - 1.f;
        sy = 2.f / float(m_sHeight) * sy - 1.f;
    }

    //moves the texture by setting its center at (scx,scy)
    void SetTcenPos(float scx, float scy){
        float tcx = float(m_tWidth) * 0.5f, tcy = float(m_tHeight) * 0.5f;
        m_tOx = scx - m_scale * (ca * tcx - sa * tcy);
        m_tOy = scy - m_scale * (sa * tcx + ca * tcy);
        if(m_angle == 0 && m_scale == 1.f){
            m_tOx = std::floor(m_tOx);
            m_tOy = std::floor(m_tOy);
        }
        compute_sO();
    }

    //moves the texture by setting its center the center of the screen
    void Tcenter(){
        SetTcenPos(m_sWidth/2, m_sHeight/2);
    }

    //move the texture area by (stx, sty)
    void Move(float stx, float sty){
        m_tOx += stx; m_tOy += sty;
        compute_sO();
    }

    //rotate the texture by ang° around origin
    void Rotate(float ang){
        m_angle += ang;
        float c=std::cos(deg2rad(ang)), s=std::sin(deg2rad(ang));
        float tox = m_tOx, toy = m_tOy;
        m_tOx = c * tox - s * toy;
        m_tOy = s * tox + c * toy;
        ca = std::cos(deg2rad(m_angle)); sa = std::sin(deg2rad(m_angle));
        compute_sO();
    }

    //rotate the texture by ang° around point(scx,scy)
    void Rotate(float ang, float scx, float scy){
        Move(-scx, -scy);
        Rotate(ang);
        Move(scx, scy);
    }

    //scale the texture by scl around origin
    void Scale(float scl){
        m_scale *= scl;
        m_tOx *= scl;
        m_tOy *= scl;
        compute_sO();
    }

    //scale the texture by scl around point(scx,scy)
    void Scale(float scl, float scx, float scy){
        Move(-scx, -scy);
        Scale(scl);
        Move(scx, scy);
    }

    void SetTscale(float scale){
        m_scale = scale;
        compute_sO();
    }

    void SetTangle(float angle){
        m_angle = angle;
        ca = std::cos(deg2rad(angle)); sa = std::sin(deg2rad(angle));
        compute_sO();
    }
private:
    void compute_sO(){
        m_sOx = -1.f/m_scale *( ca * m_tOx + sa * m_tOy);
        m_sOy = -1.f/m_scale *(-sa * m_tOx + ca * m_tOy);
    }
    float deg2rad(float deg){
        return deg * 3.14159f / 180.f;
    }
};

#endif // IMG2SCR_H
