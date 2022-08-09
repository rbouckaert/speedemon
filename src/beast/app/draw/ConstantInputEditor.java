package beast.app.draw;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.DoubleInputEditor;
import javafx.scene.control.Tooltip;
import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;

public class ConstantInputEditor extends DoubleInputEditor {
	private static final long serialVersionUID = 1L;
	
	@Override
	public Class<?> type() {
		return Function.Constant.class;
	}

	public ConstantInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
			boolean addButtons) {
		super.init(input, beastObject, itemNr, isExpandOption, addButtons);
		m_entry.setText(((Function.Constant)m_input.get()).getValue());
	}
	
	
	@Override
    protected void initEntry() {
        if (m_input.get() != null) {
            m_entry.setText(m_input.get().toString());
        }
    }

	@Override
    protected void processEntry() {
        try {
        	((Function.Constant)m_input.get()).setValue(m_entry.getText());
            validateInput();
            m_entry.requestFocus();
        } catch (Exception ex) {
            if (m_validateLabel != null) {
                m_validateLabel.setVisible(true);
                m_validateLabel.setTooltip(new Tooltip("<html><p>Parsing error: " + ex.getMessage() + ". Value was left at " + m_input.get() + ".</p></html>"));
                m_validateLabel.setColor("orange");
            }
            repaint();
        }
    }

	
}
