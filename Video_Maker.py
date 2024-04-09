import cv2, os
class main():
    def __init__(self):
        self.frame_folder = "C:\\Users\\rollo\\Documents\\GitHub\\Generative-Slicer\\Optimisation video"
        self.output_name = "video.avi"
        self.output_folder = "C:\\Users\\rollo\\Documents\\GitHub\\Generative-Slicer\\Videos"
        self.mainloop()
    def mainloop(self):
        self.files = sorted(list(os.walk(self.frame_folder))[0][2], key = lambda x: int(x.replace('.png','')[0:3]))
        print(self.files)
        print("Calibration frame address")
        print(self.frame_folder+"\\"+self.files[0])
        self.calib_frame = cv2.imread(self.frame_folder+"\\"+self.files[0])
        self.height, self.width, _  = self.calib_frame.shape

        self.video = cv2.VideoWriter(self.output_folder+"\\"+self.output_name, 0, 20, (self.width,self.height))
        for frame in self.files:
            self.video.write(cv2.imread(self.frame_folder+"\\"+frame))
        cv2.destroyAllWindows()
        self.video.release()

if __name__ == '__main__':
    main()
