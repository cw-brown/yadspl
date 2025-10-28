#include "packer_former.h"
#include "crcinterface.h"
#include <bit>

PacketFormer::PacketFormer(
        uint8_t sender_address,
        uint8_t reciever_address,
        bool are_we_source) :
        //inputBuffer(input_buffer),
        //outputBuffer(output_buffer),
        senderAddress(sender_address),
        recieverAddress(reciever_address),
        areWeSource(are_we_source) {

    // Make crc generator with polynomial for ISO 3309 in reversed format
    crcGenny = crcutil_interface::CRC::Create(0xEDB88320, 0, 32, true, 0, 0, 0, true, NULL);
    tempPacket.senderAddress = senderAddress;
    tempPacket.recieverAddress = recieverAddress;

}

bool PacketFormer::getAreWeSource() {

    return areWeSource;

}

void PacketFormer::setAreWeSource(bool are_we_source) {

    areWeSource = are_we_source;
    
}

Packet PacketFormer::formDataPacket(std::vector<uint8_t> * data, uint8_t sequence_number, uint8_t data_length) {

    crcutil_interface::UINT64 tempCRC;

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    
    tempPacket.controlCode = 0b0000;

    tempPacket.dataLength = data_length;

    tempPacket.sequenceNumber = sequence_number;

    // Form data section
    for(int i = data_length - 1; i >= 0; i--) {

        (tempPacket.data).push_back(*data->end());
        data->pop_back();

    }

    // Generate CRC and store to ERC field
    crcGenny->Compute(tempPacket.data.data(), tempPacket.dataLength, &tempCRC, NULL);
    tempPacket.erc = (uint32_t) tempCRC;

    return tempPacket;

}

Packet PacketFormer::formRetransmitPacket(std::vector<uint8_t> * data, uint8_t sequence_num, uint8_t data_length) {
    
    // Reset packet contents
    (tempPacket.data).clear();
    
    // Form initial packet with class method
    tempPacket = formDataPacket(data, sequence_num, data_length);

    // Set packet field(s)

    tempPacket.controlCode = 0b0001;

    return tempPacket;

}

Packet PacketFormer::formRetransmitPacket(Packet data_packet) {

    // Set packet field(s)
    data_packet.controlCode = 0b0001;

    return data_packet;

}

Packet PacketFormer::formBusyStartPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0010;

    tempPacket.dataLength = 0;

    tempPacket.sequenceNumber = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formBusyEndPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0011;

    tempPacket.dataLength = 0;

    tempPacket.sequenceNumber = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formCenterFrequencyPacket(const float center_freq) {

    // To hold float temporarily
    uint32_t floatConv;
    std::vector<uint8_t> floatConvBroken = {0, 0, 0, 0};

    // Reset packet contents
    (tempPacket.data).clear();

    // Cast float
    floatConv = std::bit_cast<uint32_t>(center_freq);

    // Break float into bytes
    // (floatConv >> 8*N) & 0b11111111 is a way to extract the Nth 8-bit section of floatConv from the right
    floatConvBroken[0] = ((floatConv >> 8*3) & 0b11111111);
    floatConvBroken[0] = ((floatConv >> 8*2) & 0b11111111);
    floatConvBroken[0] = ((floatConv >> 8*1) & 0b11111111);
    floatConvBroken[0] = ((floatConv >> 8*0) & 0b11111111);

    tempPacket = formDataPacket(&floatConvBroken, 0, 4);

    // Set packet field(s)
    tempPacket.controlCode = 0b0100;

    return tempPacket;

}

Packet PacketFormer::formModulationChangePacket(const uint8_t modulationType) {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0101;

    tempPacket.dataLength = 1;

    tempPacket.sequenceNumber = 0;\

    (tempPacket.data).push_back(modulationType);

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formTransactionStartPacket(const float center_freq, const uint8_t modulation_type) {
    
    // Reset packet contents
    (tempPacket.data).clear();

    // Add frequency information
    tempPacket = formCenterFrequencyPacket(center_freq);

    // Add modulation method information
    (tempPacket.data).push_back(modulation_type);

    // Set packet field(s)
    tempPacket.controlCode = 0b1000;

    tempPacket.dataLength = 5;

    tempPacket.sequenceNumber = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formTransactionRestartPacket(const float center_freq, const uint8_t modulation_type) {

    // Reset packet contents
    (tempPacket.data).clear();

    // Add transaction information
    tempPacket = formTransactionStartPacket(center_freq, modulation_type);

    // Set packet field(s)
    tempPacket.controlCode = 0b1001;

    return tempPacket;

}

Packet PacketFormer::formTransactionTransferPacket(const float center_freq, const uint8_t modulation_type) {

    // Reset packet contents
    (tempPacket.data).clear();

    // Add transaction information
    tempPacket = formTransactionStartPacket(center_freq, modulation_type);

    // Set packet field(s)
    tempPacket.controlCode = 0b1010;

    return tempPacket;

}

Packet PacketFormer::formTxTransactionDroppedPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b1011;

    tempPacket.dataLength = 0;

    tempPacket.sequenceNumber = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formTransactionEndPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b1111;

    tempPacket.dataLength = 0;

    tempPacket.sequenceNumber = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

// Reciever Packets

Packet PacketFormer::formAcknowledgePacket(uint16_t error_count, uint8_t sequence_number) {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0000;

    tempPacket.sequenceNumber = sequence_number;

    tempPacket.erc = 0;

    // Error information
    tempPacket.dataLength = 2;
    (tempPacket.data).push_back((error_count >> 8) & 0b11111111);
    (tempPacket.data).push_back((error_count >> 0) & 0b11111111);

    return tempPacket;

}

Packet PacketFormer::formRepeatPacket(uint8_t sequence_number) {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0001;

    tempPacket.sequenceNumber = sequence_number;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formBusyAcceptPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0010;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formBusyEndedPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0011;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formGeneralChannelAcknowledge() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0100;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formTransactionStartAcknowledge() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b1000;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formTransferAccept() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b1001;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formTransferDecline() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b1010;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formRxTransactionDroppedPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b1011;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}

Packet PacketFormer::formTransactionEndAcknowledge() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b1111;

    tempPacket.sequenceNumber = 0;

    tempPacket.dataLength = 0;

    tempPacket.erc = 0;

    return tempPacket;

}