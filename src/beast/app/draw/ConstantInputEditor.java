package beast.app.draw;

import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.Function;
import beast.core.Input;

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

    protected void processEntry() {
        try {
        	((Function.Constant)m_input.get()).setValue(m_entry.getText());
            m_entry.requestFocusInWindow();
        } catch (Exception ex) {
            repaint();
        }
    }

	
}
