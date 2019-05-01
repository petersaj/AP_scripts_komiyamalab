import java.io.*;
import java.awt.image.*;
import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.process.*;
import ij.plugin.PlugIn;
import ij.measure.*;

public class IJStackAlign_AP implements PlugIn
{
        public IJStackAlign_AP() {;}
        
        public ImagePlus doAlign(String cmdstr, ImagePlus source, ImagePlus target) { 
            StackReg_ tr = new StackReg_();
            int nSlices = source.getNSlices();
            for (int i = 0; i < nSlices; i++) {
                source.setSlice((i+1));
                ImagePlus sourceSlice = new ImagePlus("",source.getProcessor());
                tr.run(cmdstr, sourceSlice, target);
                ImagePlus alignedSlice = tr.getTransformedImage();
                alignedSlice.setSlice(1);
                source.setProcessor("", alignedSlice.getProcessor());
			}
            source.setSlice(1);
            return source;
        }
        
        public void run(String arg) {;}
}