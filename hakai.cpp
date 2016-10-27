#include <chrono>
#include <vector>
#include <map>
#include <forward_list>
#include <random>
#include <memory>
#include <locale>
#include <codecvt>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <DxLib.h>
#include "bmp.hpp"

double rnd(){
    static std::mt19937 mt(0x23242526);
    return (mt() % (1 << 30)) / static_cast<double>(1 << 30);
}

using namespace Eigen;
using vector2d = Eigen::Matrix<double,2,1,Eigen::DontAlign>;

// ウィンドウサイズ
const int window_width = 800;
const int window_height = 600;

const double pi = 4.0 * std::atan(1.0);

// マスク
std::unique_ptr<unsigned char, std::default_delete<unsigned char[]>> mask_data(new unsigned char[window_width * window_height]);

void clear_mask(){
    for(int i = 0; i < window_width; ++i){
        for(int j = 0; j < window_height; ++j){
            mask_data.get()[j * window_width + i] = 0xFF;
        }
    }
}

void set_mask(int x, int y){
    mask_data.get()[y * window_width + x] = 0x00;
}

//-------- FPSマネージャー
class fps_manager{
public:
    // 目標FPS
    const int target_fps = 60;

    // 修正秒率
    const int fix_time_ratio = 10;

    // 修正フレーム単位
    const int fix_frames_num = target_fps / fix_time_ratio;

    // スリープタイム
    const float sleep_time = 1000.0f / target_fps;

    // 表示
    void draw_sleep_time() const{
        char str[7] = { 0 };
        str[0] = '0' + (static_cast<int>(fps.real_sleep_time / 10) % 10);
        str[1] = '0' + (static_cast<int>(fps.real_sleep_time) % 10);
        str[2] = '.';
        str[3] = '0' + (static_cast<int>(fps.real_sleep_time * 10) % 10);
        str[4] = '0' + (static_cast<int>(fps.real_sleep_time * 100) % 10);
        str[5] = '0' + (static_cast<int>(fps.real_sleep_time * 1000) % 10);
        DrawString(0, 0, str, GetColor(0x00, 0x00, 0x00));
    }

    // 処理
    void operator ()(){
        ++frame_count;
        Sleep(static_cast<int>(real_sleep_time));
        if(frame_count >= fix_frames_num){
            std::chrono::milliseconds real = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
            real_sleep_time = sleep_time * (sleep_time * fix_frames_num) / static_cast<float>(real.count());
            if(real_sleep_time >= sleep_time){
                real_sleep_time = sleep_time;
            }
            start = std::chrono::high_resolution_clock::now();
            frame_count = 0;
        }
    }

    float real_sleep_time = 16;

private:
    int frame_count = 0;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
} fps;

//-------- ビット配列
class bitarray{
public:
    bitarray() = default;
    bitarray(std::size_t n) :
        size_(n / (sizeof(unsigned int) * 8) + ((n % (sizeof(unsigned int) * 8)) > 0 ? 1 : 0)),
        data(size_)
    {}

    bool get(std::size_t i) const{
        return ((data[i / (sizeof(unsigned int) * 8)] >> (i % (sizeof(unsigned int) * 8))) & 1) == 1;
    }

    void set(std::size_t i, bool x){
        unsigned int &v = data[i / (sizeof(unsigned int) * 8)];
        v = x ? v | (1 << (i % (sizeof(unsigned int) * 8))) : v & ~(1 << (i % (sizeof(unsigned int) * 8)));
    }

    void resize(std::size_t n){
        size_ = n / (sizeof(unsigned int) * 8) + (n % (sizeof(unsigned int) * 8) > 0 ? 1 : 0);
        data.resize(size_);
    }

    void clear(){
        size_ = 0;
        data.clear();
    }

    void zero_init(){
        for(auto &i : data){
            i = 0;
        }
    }

    std::size_t size() const{
        return size_;
    }

private:
    std::vector<unsigned int> data;
    std::size_t size_;
};

//-------- テクスチャ
class texture{
public:
    void make(int w, int h){
        handle = MakeScreen(w, h, TRUE);
        SetDrawScreen(handle);
        SetTransColor(0x00, 0xFF, 0xFF);
        DrawBox(0, 0, w, h, GetColor(0xFF, 0xFF, 0xFF), TRUE);
    }

    int get() const{
        return handle;
    }

    int width() const{
        return width_;
    }

    int height() const{
        return height_;
    }

private:
    int handle;
    int width_ = 0, height_ = 0;
};

//-------- 影
class collision{
public:
    void make(int width, int height){
        width_ = width;
        height_ = height;
        n_width_ = (std::max)(width, height) * 3 / 2;
        n_height_ = (std::max)(width, height) * 3 / 2;
        ba.resize(width * height);
        rotated.resize(n_width_ * n_height_);
        ba.zero_init();
        rotated.zero_init();
    }

    bool get(int x, int y) const{
        return ba.get(y * width_ + x);
    }

    void set(int x, int y, bool a){
        ba.set(y * width_ + x, a);
    }

    void rotete(double rad){
        rotated.zero_init();
        vector2d center, n_center;
        center << width_ / 2, height_ / 2;
        n_center << n_width_ / 2, n_height_ / 2;
        for(int ny = 0; ny < n_height_; ++ny){
            for(int nx = 0; nx < n_width_; ++nx){
                vector2d p;
                p << nx, ny;
                p -= n_center;
                double theta = std::atan2(p[1], p[0]);
                double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
                vector2d q;
                q << std::cos(theta - rad) * r, std::sin(theta - rad) * r;
                q += center;
                if(
                    q[0] >= 0 &&
                    q[0] < width_ &&
                    q[1] >= 0 &&
                    q[1] < height_ &&
                    get(static_cast<int>(q[0]), static_cast<int>(q[1]))
                ){
                    rotated.set(ny * n_width_ + nx, true);
                }
            }
        }
    }

    int width() const{
        return width_;
    }

    int height() const{
        return height_;
    }

    int n_width() const{
        return n_width_;
    }

    int n_height() const{
        return n_height_;
    }

    bool collision_get(int x, int y){
        return rotated.get(y * n_width_ + x);
    }

private:
    bitarray ba;
    bitarray rotated;
    int width_ = 0, height_ = 0;
    int n_width_ = 0, n_height_ = 0;
};

//-------- 質点
class mass_point{
public:
    static std::vector<mass_point> point_list;

    mass_point(){
        coord << 0.0, 0.0;
        speed << 0.0, 0.0;
    }

    mass_point(const mass_point &) = default;
    ~mass_point() = default;

    void init_speed(){
        vector2d window_center, p;
        window_center << window_width / 2, window_height / 2;
        p = coord - window_center;
        double theta = std::atan2(p[1], p[0]) + (pi / (20.0 + (rnd() * 8.0 - 4.0))) * (rnd() < 0.5 ? +1 : -1);
        double r = 1.0 + rnd() * 0.5;
        speed <<
            std::cos(theta) * r,
            std::sin(theta) * r;
        rad_speed = (pi / 22.0) * (rnd() < 0.5 ? +1 : -1);
    }

    void update_speed(){
        rad += rad_speed;
    }

    void update_position(){
        coord = coord + speed;
    }

    // 描画する
    void draw(){
        SetDrawBlendMode(DX_BLENDMODE_MULA, 0xFF) ;
        DrawRotaGraph(
            static_cast<int>(coord[0]),
            static_cast<int>(coord[1]),
            1.0, rad, tex.get(), TRUE
        );
    }

    // 質量
    double mass = 0.0;

    // 座標
    vector2d coord;

    // 速度
    vector2d speed;

    // 半径
    double radius = 0.0;

    // 色
    unsigned int color;

    // テクスチャ
    texture tex;

    // 衝突判定用の2値影
    collision sha;

    // 角度
    double rad = 0.0;
    double rad_speed = 0.0;

    // 軌道
    std::vector<vector2d> orbit;
};

std::vector<mass_point> mass_point::point_list;

void init_texture(const wchar_t *path){
    tt_legacy::bmp bmp(path);

    auto make_texture = [&](int x, int y){
        int or_x = x, or_y = y;
        int min_x = x, min_y = y;
        int max_x = 0, max_y = 0;
        
        std::vector<std::pair<int, int>> target = { std::make_pair(x, y) }, scaned;
        
        while(!target.empty()){
            std::pair<int, int> a = target.back();
            target.pop_back();
            scaned.push_back(std::make_pair(a.first - or_x, a.second - or_y));
            
            min_x = (std::min)(a.first - or_x, min_x);
            min_y = (std::min)(a.second - or_y, min_y);

            max_x = (std::max)(a.first - or_x, max_x);
            max_y = (std::max)(a.second - or_y, max_y);
            
            bmp.clr(a.first, a.second, bmp.rgb(0xFF, 0xFF, 0xFF));
            
            if(a.first + 1 < bmp.width() && bmp.clr(a.first + 1, a.second) == bmp.rgb(0x00, 0x00, 0x00)){
                target.push_back(std::make_pair(a.first + 1, a.second));
            }
            if(a.first - 1 > 0 && bmp.clr(a.first - 1, a.second) == bmp.rgb(0x00, 0x00, 0x00)){
                target.push_back(std::make_pair(a.first - 1, a.second));
            }
            if(a.second + 1 < bmp.height() && bmp.clr(a.first, a.second + 1) == bmp.rgb(0x00, 0x00, 0x00)){
                target.push_back(std::make_pair(a.first, a.second + 1));
            }
            if(a.second - 1 > 0 && bmp.clr(a.first, a.second - 1) == bmp.rgb(0x00, 0x00, 0x00)){
                target.push_back(std::make_pair(a.first, a.second - 1));
            }
        }

        {
            or_x += min_x;
            or_y += min_y;
            max_x += -min_x;
            max_y += -min_y;
            max_x = max_x / 2 + (max_x & 1);
            max_y = max_y / 2 + (max_y & 1);
            for(auto &i : scaned){
                i.first += -min_x - max_x;
                i.second += -min_y - max_y;
            }
            min_x = -max_x;
            min_y = -max_y;
        }

        int width = max_x - min_x + 1, height = max_y - min_y + 1;
        mass_point m;
        m.mass = 1.0;
        m.coord <<
            window_width / 2 - bmp.width() / 2 + or_x + width / 2 + (width & 1),
            window_height / 2 - bmp.height() / 2 + or_y + height / 2 + (height & 1);
        m.speed << 0.0, 0.0;
        m.radius = std::sqrt(width * width + height * height) / 2.0;
        m.tex.make(width, height);
        m.sha.make(width, height);
        SetDrawScreen(m.tex.get());
        for(std::size_t i = 0; i < scaned.size(); ++i){
            int x_ = scaned[i].first;
            int y_ = scaned[i].second;
            DrawPixel(-min_x + x_, -min_y + y_, GetColor(0x00, 0x00, 0x00));
            m.sha.set(-min_x + x_, -min_y + y_, true);
        }
        mass_point::point_list.push_back(m);
    };

    for(int y = 0; y < bmp.height(); ++y){
        for(int x = 0; x < bmp.width(); ++x){
            if(bmp.clr(x, y) == bmp.rgb(0x00, 0x00, 0x00)){
                make_texture(x, y);
            }
        }
    }
}

LRESULT CALLBACK window_proc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp){
    switch(msg){
    default:
        return DefWindowProc(hwnd, msg, wp, lp);
    }

    return 0;
}

int WINAPI WinMain(HINSTANCE handle, HINSTANCE prev_handle, LPSTR lp_cmd, int n_cmd_show){
    // init DxLib
    if(
        SetOutApplicationLogValidFlag(FALSE) != 0 ||
        ChangeWindowMode(TRUE) != DX_CHANGESCREEN_OK ||
        SetGraphMode(window_width, window_height, 32) != DX_CHANGESCREEN_OK ||
        SetMainWindowText("hakai") != 0 ||
        DxLib_Init() != 0
    ){ return -1; }

    SetHookWinProc(window_proc);
    {
        namespace fs = boost::filesystem;
        fs::path path("./");
        BOOST_FOREACH(const fs::path &p, std::make_pair(fs::directory_iterator(path), fs::directory_iterator())){
            if(!fs::is_directory(p)){
                std::wstring name = p.filename().c_str();
                if(name.substr(name.size() - 4, 4) == std::wstring(L".bmp")){
                    init_texture(name.c_str());
                }
            }
        }
    }

    SetDrawScreen(DX_SCREEN_BACK);
    SetTransColor(0x00, 0xFF, 0xFF);

    bool flag = true;

    while(true){
        if(ProcessMessage() == -1){
            break;
        }

        // draw
        //SetDrawScreen(main_graphic_handle);
        SetDrawBlendMode(DX_BLENDMODE_NOBLEND, 0xFF);
        DrawBox(0, 0, window_width, window_height, GetColor(0xFF, 0xFF, 0xFF), TRUE);

        for(mass_point &i : mass_point::point_list){
            i.draw();
        }

        ScreenFlip();

        if(flag){
            WaitKey();
            flag = false;

            for(mass_point &i : mass_point::point_list){
                i.init_speed();
            }
        }

        for(mass_point &i : mass_point::point_list){
            i.update_speed();
        }

        for(mass_point &i : mass_point::point_list){
            i.update_position();
        }

        // fps
        fps();
    }

    DxLib_End();
    return 0;
}
